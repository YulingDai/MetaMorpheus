using Chemistry;
using EngineLayer.spectralLibrarySearch;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ClassicSearch
{
    public class ClassicSearchEngine : MetaMorpheusEngine
    {
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Protein> Proteins;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly List<SilacLabel> SilacLabels;
        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;
        private readonly SpectralLibrarySearchResults[] DotPeptideSpectralMatches;

        public Dictionary<String, Spectrum> SpectralLibraryDictionary { get; set; }

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications, List<Modification> fixedModifications, List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel,
            List<Protein> proteinList, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            PeptideSpectralMatches = globalPsms;
            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
            VariableModifications = variableModifications;
            FixedModifications = fixedModifications;
            SilacLabels = silacLabels;
            if (startLabel != null || endLabel != null) //else it's null
            {
                TurnoverLabels = (startLabel, endLabel);
            }
            Proteins = proteinList;
            SearchMode = searchMode;
            DotPeptideSpectralMatches = new SpectralLibrarySearchResults[PeptideSpectralMatches.Length];
        }



        protected override MetaMorpheusEngineResults RunSpecific()
        {
            int xyyyy = 0;
            Status("Getting ms2 scans...");

            double proteinsSearched = 0;
            int oldPercentProgress = 0;

            // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
            var myLocks = new object[PeptideSpectralMatches.Length];
            for (int i = 0; i < myLocks.Length; i++)
            {
                myLocks[i] = new object();
            }

            Status("Performing classic search...");

            if (Proteins.Any())
            {
                int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
                int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
                Parallel.ForEach(threads, (i) =>
                {
                    var peptideTheorProducts = new List<Product>();

                    for (; i < Proteins.Count; i += maxThreadsPerFile)
                    {
                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { return; }

                        // digest each protein into peptides and search for each peptide in all spectra within precursor mass tolerance
                        foreach (PeptideWithSetModifications peptide in Proteins[i].Digest(CommonParameters.DigestionParams, FixedModifications, VariableModifications, SilacLabels, TurnoverLabels))
                        {

                            peptide.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);

                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(peptide.MonoisotopicMass, SearchMode))

                            {

                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters);

                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);
                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;
                                xyyyy++;

                                // this is thread-safe because even if the score improves from another thread writing to this PSM,
                                // the lock combined with AddOrReplace method will ensure thread safety
                                if (meetsScoreCutoff)
                                {

                                    // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                    lock (myLocks[scan.ScanIndex])
                                    {
                                        bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (thisScore - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                                        if (scoreImprovement)
                                        {
                                            if (PeptideSpectralMatches[scan.ScanIndex] == null)
                                            {
                                                PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(peptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, CommonParameters, matchedIons, 0);
                                            }
                                            else
                                            {
                                                PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(peptide, thisScore, scan.Notch, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                                            }

                                            //if (SpectralLibraryDictionary != null && SpectralLibraryDictionary.ContainsKey(peptide.FullSequence))
                                            //{
                                            //    var spectrumMz = SpectralLibraryDictionary[peptide.FullSequence].Peaks.Select(p => p.NeutralMass).ToArray();
                                            //    var spectrumIntensity = SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity.Select(p => p.Value).ToArray();
                                            //    var librarySpectrum = new MzSpectrum(spectrumMz, spectrumIntensity, false);
                                            //    List<MatchedFragmentIon> spectrumMatchedIons = MatchSpectrumPeaks(scan.TheScan, SpectralLibraryDictionary[peptide.FullSequence].Peaks, CommonParameters);
                                            //    var matchedMz = spectrumMatchedIons.Select(p => p.Mz).ToArray();
                                            //    var matchedIntensity = spectrumMatchedIons.Select(p => p.Intensity).ToArray();
                                            //    var matchedpectrum = new MzSpectrum(matchedMz, matchedIntensity, false);
                                            //    double spectrumDotProductScore = matchedpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                            //    //double spectrumDotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                            //    PeptideSpectralMatches[scan.ScanIndex].Score = thisScore;
                                            //    PeptideSpectralMatches[scan.ScanIndex].ThisSpectrumMatchScore = spectrumDotProductScore;
                                            //}
                                        }
                                    }
                                }

                                if (SpectralLibraryDictionary != null && SpectralLibraryDictionary.ContainsKey(peptide.FullSequence))
                                {
                                    var spectrumMz = SpectralLibraryDictionary[peptide.FullSequence].Peaks.Select(p => p.NeutralMass).ToArray();
                                    var spectrumIntensity = SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity.Select(p => p.Value).ToArray();
                                    var librarySpectrum = new MzSpectrum(spectrumMz, spectrumIntensity, false);
                                    List<MatchedFragmentIon> spectrumMatchedIons = MatchSpectrumPeaks(scan.TheScan, SpectralLibraryDictionary[peptide.FullSequence].Peaks, CommonParameters);
                                    var matchedMz = spectrumMatchedIons.Select(p => p.Mz).ToArray();
                                    var matchedIntensity = spectrumMatchedIons.Select(p => p.Intensity).ToArray();
                                    var matchedpectrum = new MzSpectrum(matchedMz, matchedIntensity, false);
                                    //double spectrumDotProductScore = matchedpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                    double spectrumDotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                    //if (peptide.FullSequence == "RPYESSR")
                                    //{
                                        double spectrumSharedDotProductScore = this.CalculateSharedDotProductSimilarity(scan.TheScan.TheScan.MassSpectrum, librarySpectrum, CommonParameters.ProductMassTolerance);
                                    //}
                                    //double spectrumSharedDotProductScore = this.CalculateSharedDotProductSimilarity(scan.TheScan.TheScan.MassSpectrum, librarySpectrum, CommonParameters.ProductMassTolerance);
                                    if(thisScore < 5)
                                    {
                                        Console.WriteLine(thisScore + "  " + spectrumDotProductScore + "   " + spectrumSharedDotProductScore);
                                    }
                                    
                                    //Console.WriteLine(peptide.FullSequence +"   " +thisScore + "  " + spectrumDotProductScore + "   " );
                                    //if (DotPeptideSpectralMatches[scan.ScanIndex] == null || spectrumDotProductScore > DotPeptideSpectralMatches[scan.ScanIndex].Score)
                                    //{
                                    //    DotPeptideSpectralMatches[scan.ScanIndex] = new SpectralLibrarySearchResults(peptide.FullSequence, spectrumDotProductScore, scan.ScanIndex);

                                    //}


                                    //if (!meetsScoreCutoff && spectrumDotProductScore>0.8)
                                    //{
                                    //    Console.WriteLine(thisScore + "     " + spectrumDotProductScore);
                                    //}
                                }




                                //if (SpectralLibraryDictionary != null && SpectralLibraryDictionary.ContainsKey(peptide.FullSequence))
                                //{

                                //      List<MatchedFragmentIon> spectrumMatchedIons = MatchSpectrumPeaks(scan.TheScan, SpectralLibraryDictionary[peptide.FullSequence].Peaks, CommonParameters);
                                //     double thisSpectrumMatchScore = CalculatePeptideScore(scan.TheScan.TheScan, spectrumMatchedIons);

                                //    var spectrumMz = SpectralLibraryDictionary[peptide.FullSequence].Peaks.Select(p => p.NeutralMass).ToArray();
                                //    var spectrumIntensity = SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity.Select(p => p.Value).ToArray();
                                //    var librarySpectrum = new MzSpectrum(spectrumMz, spectrumIntensity, false);
                                //    double spectrumDotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                //    PeptideSpectralMatches[scan.ScanIndex].Score = thisScore;
                                //    PeptideSpectralMatches[scan.ScanIndex].ThisSpectrumMatchScore = spectrumDotProductScore;

                                //Console.WriteLine(scan.ScanIndex +"   " +  peptide.FullSequence + "   " + thisSpectrumMatchScore + "   " + thisScore + "   " + spectrumDotProductScore);
                                //var theoreticalMz = SpectralLibraryDictionary[peptide.FullSequence].Peaks.Select(p => p.NeutralMass).ToArray();
                                //var spectrumIntensity = SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity.Select(p => p.Value).ToArray();
                                //var librarySpectrum = new MzSpectrum(spectrumMz, spectrumIntensity, false);
                                //double spectrumDotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                //var mz = SpectralLibraryDictionary[peptide.FullSequence].Peaks.Select(p => p.NeutralMass).ToArray();
                                //var intensity = SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity.Select(p => p.Value).ToArray();
                                //if (peptide.FullSequence == "SKDVTDSATTKK")
                                //{
                                //    Console.WriteLine(scan.ScanIndex);
                                //    Console.WriteLine(peptide.FullSequence + " " + thisSpectrumMatchScore + "  " + spectrumMatchedIons.Count + "  " + SpectralLibraryDictionary[peptide.FullSequence].Peaks.Count + "  " + matchedIons.Count);
                                //    Console.WriteLine(spectrumMatchedIons.Count);
                                //    foreach (var x in spectrumMatchedIons)
                                //    {
                                //        Console.WriteLine(x.Intensity);
                                //    }
                                //    Console.WriteLine(matchedIons.Count);
                                //    foreach (var x in matchedIons)
                                //    {
                                //        Console.WriteLine(x.Intensity);
                                //    }
                                //    //Console.WriteLine(peptideTheorProducts.Count);
                                //    //foreach (var x in peptideTheorProducts)
                                //    //{
                                //    //    Console.WriteLine(x.NeutralMass + "\t" + x.ProductType + x.FragmentNumber);
                                //    //}
                                //    Console.WriteLine(SpectralLibraryDictionary[peptide.FullSequence].Peaks.Count);
                                //    foreach (var x in SpectralLibraryDictionary[peptide.FullSequence].PeaksWithIntensity)
                                //    {
                                //        Console.Write( x.Value+ "\t" + x.Key);
                                //    }
                                //    console.WriteLine(scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray().Length);
                                //    foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.Peaks).ToArray())
                                //    {
                                //        Console.WriteLine( x.First());

                                //    }
                                //    //foreach (var x in scan.TheScan.TheScan.MassSpectrum.YArray)
                                //    //{                                        
                                //    Console.WriteLine(x + "\t");
                                //    //}
                                //    //foreach (var x in scan.TheScan.TheScan.MassSpectrum.XArray)
                                //    //{
                                //    //    Console.WriteLine(x + "\t");
                                //    //}
                                //    //    //Console.WriteLine(peptide.FullSequence + scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray());
                                //    //    //scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
                                //    //Console.WriteLine(scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray().Length);
                                //    //foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass.ToMz(p.Charge)).ToArray())
                                //    //{
                                //    //    Console.WriteLine(x + "\t");
                                //    //}
                                //}
                                //    Console.WriteLine(peptide.FullSequence);

                                //Console.WriteLine(scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray().Length);
                                //scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
                                //foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray())
                                //{
                                //    Console.WriteLine(x);
                                //}
                                //foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass.ToMz(p.Charge)).ToArray())
                                //{
                                //    Console.WriteLine(x);
                                //}

                                //Console.WriteLine(intensity.Length);
                                //intensity.ToList().ForEach(b => Console.WriteLine(b + "\t"));
                                //Console.WriteLine(mz.Length);
                                //mz.ToList().ForEach(b => Console.WriteLine(b + "\t"));


                                //var librarySpectrum = new MzSpectrum(mz, intensity, false);
                                //double spectrumDotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(librarySpectrum, CommonParameters.ProductMassTolerance);
                                //Console.WriteLine(peptide.FullSequence + "     " + scan.ScanIndex + "     " +  spectrumDotProductScore + "     " + thisScore + "      " + thisSpectrumMatchScore);
                                //}



                                //if (!meetsScoreCutoff && spectrumDotProductScore > 0.2)
                                //{
                                //    Console.Write(spectrumDotProductScore + "\t");
                                //}
                                //Console.WriteLine(spectrumDotProductScore + "\t" + thisScore + "\t" + thisSpectrumMatchScore);
                                //if (spectrumMatchedIons.Count != matchedIons.Count)
                                // { 
                                //Console.WriteLine(peptide.FullSequence + " " + thisSpectrumMatchScore + "  " + spectrumMatchedIons.Count + "  " + SpectralLibraryDictionary[peptide.FullSequence].Peaks.Count + "  " + matchedIons.Count);
                                //    Console.WriteLine(peptide.FullSequence);
                                //    //Console.WriteLine(peptide.FullSequence + " " + thisSpectrumMatchScore +" " + thisScore +  "  " + spectrumMatchedIons.Count + "  " + SpectralLibraryDictionary[peptide.FullSequence].Peaks.Count + "  " + matchedIons.Count);
                                //    foreach (var x in spectrumMatchedIons)
                                //    {

                                //        Console.Write(x.Mz + "\t");
                                //    }
                                //    foreach (var x in matchedIons)
                                //    {
                                //        Console.Write(x.Mz + "\t");
                                //    }
                                //    foreach (var x in peptideTheorProducts)
                                //    {
                                //        Console.Write(x.NeutralMass + "\t" + x.ProductType + x.FragmentNumber);
                                //    }
                                //    foreach (var x in SpectralLibraryDictionary[peptide.FullSequence].Peaks)
                                //    {
                                //        Console.Write(x.NeutralMass + "\t");
                                //    }
                                //    foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray())
                                //    {
                                //        Console.Write(x + "\t");
                                //    }
                                //    //Console.WriteLine(peptide.FullSequence + scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray());
                                //    //scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
                                //    foreach (var x in scan.TheScan.ExperimentalFragments.Select(p => p.MonoisotopicMass.ToMz(p.Charge)).ToArray())
                                //    {
                                //        Console.Write(x + "\t");
                                //    }
                                //    Console.WriteLine(peptide.FullSequence);
                                //}
                                //Console.WriteLine( peptide.FullSequence + " " + thisSpectrumMatchScore + "  " +spectrumMatchedIons.Count + "  " + SpectralLibraryDictionary[peptide.FullSequence].Peaks.Count + "  " + matchedIons.Count);
                                //bool meetSpectrumScoreCutoff = thisSpectrumMatchScore >= CommonParameters.ScoreCutoff;
                                //    if (meetSpectrumScoreCutoff)
                                //    {
                                //        // valid hit (met the cutoff score); lock the scan to prevent other threads from accessing it
                                //        lock (myLocks[scan.ScanIndex])
                                //        {
                                //            bool scoreImprovement = PeptideSpectralMatches[scan.ScanIndex] == null || (thisSpectrumMatchScore - PeptideSpectralMatches[scan.ScanIndex].RunnerUpScore) > -PeptideSpectralMatch.ToleranceForScoreDifferentiation;

                                //            if (scoreImprovement)
                                //            {
                                //                if (PeptideSpectralMatches[scan.ScanIndex] == null)
                                //                {
                                //                    PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(peptide, scan.Notch, thisSpectrumMatchScore, scan.ScanIndex, scan.TheScan, CommonParameters, spectrumMatchedIons, 0);
                                //                    PeptideSpectralMatches[scan.ScanIndex].Score = thisScore;
                                //                    PeptideSpectralMatches[scan.ScanIndex].ThisSpectrumMatchScore = thisSpectrumMatchScore;
                                //                }
                                //                else
                                //                {
                                //                    PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(peptide, thisSpectrumMatchScore, scan.Notch, CommonParameters.ReportAllAmbiguity, spectrumMatchedIons, 0);
                                //                    PeptideSpectralMatches[scan.ScanIndex].Score = thisScore;
                                //                    PeptideSpectralMatches[scan.ScanIndex].ThisSpectrumMatchScore = thisSpectrumMatchScore;
                                //            }
                                //            }
                                //        }
                                //    }

                                //double dotProductScore = scan.TheScan.TheScan.MassSpectrum.CalculateDotProductSimilarity(, CommonParameters.ProductMassTolerance);

                                //}

                            }
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        proteinsSearched++;
                        var percentProgress = (int)((proteinsSearched / Proteins.Count) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", NestedIds));
                        }
                    }
                });
            }
            // Console.WriteLine("peptideMatch");
            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
                //double intensitySum = psm.MatchedFragmentIons.Select(m => m.Intensity).Sum();
                //Console.WriteLine(psm.FullSequence + "  " + psm.Score + "  " + psm.ScanIndex);
                //foreach(var x in psm.MatchedFragmentIons)
                //{
                //    Console.WriteLine(x.Mz + "  " + x.Intensity/ intensitySum + "  " + x.NeutralTheoreticalProduct.ProductType.ToString() + x.NeutralTheoreticalProduct.FragmentNumber);
                //}

            }
            //Console.WriteLine("dotMatch");
            foreach (var x in DotPeptideSpectralMatches.Where(p => p != null))
            {
                Console.WriteLine(x.Peptide + "  " + x.Score + "  " + x.Index);
            }
            return new MetaMorpheusEngineResults(this);

        }
        public double CalculateSharedDotProductSimilarity(MzSpectrum thisSpectrum, MzSpectrum spectrumToCompare, Tolerance tolerance)
        {
            //get arrays of m/zs and intensities
            double[] mz1 = thisSpectrum.XArray;
            double[] intensity1 = thisSpectrum.YArray;

            double[] mz2 = spectrumToCompare.XArray;
            double[] intensity2 = spectrumToCompare.YArray;

            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();
            int i = 0; //iterate through mz1
            int j = 0; //iterate through mz2

            //find where peaks match
            while (i != mz1.Length && j != mz2.Length)
            {
                double one = mz1[i];
                double two = mz2[j];
                if (tolerance.Within(one, two)) //if match
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(intensity2[j]);
                    i++;
                    j++;
                }
                else if (one > two)
                {
                    //vector1.Add(0);
                    //vector2.Add(intensity2[j]);
                    j++;
                }
                else //two>one
                {
                    //vector1.Add(intensity1[i]);
                    //vector2.Add(0);
                    i++;
                }
            }
            //wrap up leftover peaks
            //for (; i < mz1.Length; i++)
            //{
            //    vector1.Add(intensity1[i]);
            //    vector2.Add(0);
            //}
            //for (; j < mz2.Length; j++)
            //{
            //    vector1.Add(0);
            //    vector2.Add(intensity2[j]);
            //}

            //numerator of dot product
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator of dot product
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            //return dot product
            return numerator / denominator;
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, MassDiffAcceptor searchMode)
        {
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum);
                if (scanIndex < ArrayOfSortedMS2Scans.Length)
                {
                    var scanMass = MyScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var scan = ArrayOfSortedMS2Scans[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(scan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == ArrayOfSortedMS2Scans.Length)
                        {
                            break;
                        }

                        scanMass = MyScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum)
        {
            int index = Array.BinarySearch(MyScanPrecursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }


    }
}