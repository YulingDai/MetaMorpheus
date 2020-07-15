using EngineLayer.spectralLibrarySearch;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

//namespace EngineLayer.ClassicSearch
//{
//    public class ClassicSearchForLibraryEngine : MetaMorpheusEngine
//    {
//        private readonly MassDiffAcceptor SearchMode;
//        private readonly List<Protein> Proteins;
//        private readonly List<Modification> FixedModifications;
//        private readonly List<Modification> VariableModifications;
//        private readonly List<SilacLabel> SilacLabels;
//        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
//        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
//        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
//        private readonly double[] MyScanPrecursorMasses;
//        private readonly double?[] SpectralLibraryPrecursorMasses;

//        public Spectrum[] Sorted_spectralLibrary { get; set; }
//        public SpectralLibrarySearchResults[] SpectralLibrarySearchResults { get; set; }
      
//        public new CommonParameters CommonParameters { get; set; }

//        public ClassicSearchForLibraryEngine(SpectralLibrarySearchResults[] spectralLibrarySearchResults, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
//            List<Modification> variableModifications, List<Modification> fixedModifications, List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel,
//           Spectrum[] spectralLibrary, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds)
//            : base(commonParameters, fileSpecificParameters, nestedIds)
//        {
//            SpectralLibrarySearchResults = spectralLibrarySearchResults;
//            ArrayOfSortedMS2Scans = arrayOfSortedMS2Scans;
//            MyScanPrecursorMasses = arrayOfSortedMS2Scans.Select(b => b.PrecursorMass).ToArray();
//            VariableModifications = variableModifications;
//            FixedModifications = fixedModifications;
//            SilacLabels = silacLabels;

//            Sorted_spectralLibrary = SortedByMzSpectralLibrary(spectralLibrary);

//            SpectralLibraryPrecursorMasses = Sorted_spectralLibrary.Select(b => b.precursorMz).ToArray();

//            if (startLabel != null || endLabel != null) //else it's null
//            {
//                TurnoverLabels = (startLabel, endLabel);
//            }
//            Sorted_spectralLibrary = Sorted_spectralLibrary;
//            SearchMode = searchMode;
//        }

        //protected override MetaMorpheusEngineResults RunSpecific()
        //{
        //    Status("Getting ms2 scans...");

        //    double scanSearched = 0;
        //    int oldPercentProgress = 0;

        //    // one lock for each MS2 scan; a scan can only be accessed by one thread at a time
        //    var myLocks = new object[SpectralLibrarySearchResults.Length];
        //    for (int i = 0; i < myLocks.Length; i++)
        //    {
        //        myLocks[i] = new object();
        //    }


        //    Status("Performing classic search...");

        //    int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
        //    int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
        //    Parallel.ForEach(threads, (i) =>
        //    {
        //        var peptideTheorProducts = new List<Product>();

        //        for (; i < Proteins.Count; i += maxThreadsPerFile)
        //        {
        //            // Stop loop if canceled
        //            if (GlobalVariables.StopLoops) { return; }

        //            for (int i = 0; i < ArrayOfSortedMS2Scans.Length; i++)
        //            {
        //                var expecrimentalScan = ArrayOfSortedMS2Scans[i];
        //                //SpectralLibrarySearchResults[i] = new SpectralLibrarySearchResults(expecrimentalSpectrum, i, new List<SpectralLibrarayMatch>());

        //                foreach (Spectrum matchedSpectrum in GetAcceptableLibrarySpectrums(Convert.ToDouble(expecrimentalScan.PrecursorMonoisotopicPeakMz), SearchMode))
        //                {
        //                    //List<PeaksInformationFromSpectrum> matchedPeaks = MatchSpectrumPeaks(expecrimentalSpectrum, matchedSpectrum);

        //                    double thisScore = CalculateSpectrumScore(matchedSpectrum, matchedPeaks);

        //                    bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

        //                    if (meetsScoreCutoff)
        //                    {

        //                        SpectralLibrarySearchResults[i].SpectralLibrarayMatchs.Add(new SpectralLibrarayMatch(matchedSpectrum, thisScore, matchedPeaks));

        //                    }

        //                }

        //            }
                
        //            // report search progress (proteins searched so far out of total proteins in database)
        //            scanSearched++;
        //            var percentProgress = (int)((scanSearched / ArrayOfSortedMS2Scans.Length) * 100);

        //            if (percentProgress > oldPercentProgress)
        //            {
        //                oldPercentProgress = percentProgress;
        //                ReportProgress(new ProgressEventArgs(percentProgress, "Performing classic search... ", NestedIds));
        //            }
        //        }
        //    }

        //    return new MetaMorpheusEngineResults(this);
        //}

        //helper method for sorting
//        static int Compare(Spectrum x, Spectrum y)
//        {
//            if (x.precursorMz < y.precursorMz)
//            {
//                return 1;
//            }
//            else if (x.precursorMz > y.precursorMz)
//            {
//                return -1;
//            }
//            else
//            {
//                return 0;
//            }
//        }


//        public Spectrum[] SortedByMzSpectralLibrary(Spectrum[] spectralLibrary)
//        {
//            Spectrum temp;
//            for (int i = 0; i < spectralLibrary.Length; i++)
//            {
//                for (int j = i + 1; j < spectralLibrary.Length; j++)
//                {
//                    if (Compare(spectralLibrary[i], spectralLibrary[j]) < 0)
//                    {
//                        temp = spectralLibrary[i];
//                        spectralLibrary[i] = spectralLibrary[j];
//                        spectralLibrary[j] = temp;
//                    }
//                }
//            }
//            return spectralLibrary;
//        }

//        public List<PeaksInformationFromSpectrum> MatchSpectrumPeaks(Spectrum experimentalSpectrum, Spectrum matchedSpectrumInLibrary, CommonParameters commonParameters)
//        {
//            var matchedSpectrumPeaks = new List<PeaksInformationFromSpectrum>();

//            for (int i = 0; i < matchedSpectrumInLibrary.Peaks.Length; i++)
//            {
//                var matchedPeak = matchedSpectrumInLibrary.Peaks[i];
//                for (int j = 0; j < experimentalSpectrum.Peaks.Length; j++)
//                {
//                    var experimentalPeak = experimentalSpectrum.Peaks[j];
//                    if (Math.Abs(experimentalPeak.Mz - matchedPeak.Mz) <= FragmentTolerance)
//                    {
//                        var matchedSpectrumPeak = new PeaksInformationFromSpectrum(matchedPeak.Mz, matchedPeak.Intensity);

//                        matchedSpectrumPeaks.Add(matchedSpectrumPeak);
//                    }
//                }

//            }
//            return matchedSpectrumPeaks;
//        }


//        public static double CalculateSpectrumScore(Spectrum thisSpectrum, List<PeaksInformationFromSpectrum> matchedPeaks)
//        {
//            double score = 0;
//            for (int i = 0; i < matchedPeaks.Count; i++)
//            {
//                score += 1 + matchedPeaks[i].Intensity;
//            }
//            return score;
//        }

//        private List<Spectrum> GetAcceptableLibrarySpectrums(double experimentalPreMz, MassDiffAcceptor searchMode)
//        {
//            var acceptableLibrarySpectrums = new List<Spectrum>();
//            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromObservedMass(experimentalPreMz).ToList())
//            {
//                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;

//                int index = GetFirstScanWithMassOverOrEqualInSpectralLibrary(allowedInterval.Minimum);
//                if (index < SpectralLibraryPrecursorMasses.Length)
//                {
//                    var precursorMass = SpectralLibraryPrecursorMasses[index];
//                    while (precursorMass <= allowedInterval.Maximum)
//                    {
//                        Spectrum matchedSpectrum = Sorted_spectralLibrary[index];
//                        acceptableLibrarySpectrums.Add(matchedSpectrum);
//                        index++;
//                        if (index == Sorted_spectralLibrary.Length)
//                        {
//                            break;
//                        }

//                        precursorMass = SpectralLibraryPrecursorMasses[index];
//                    }
//                }
//            }
//            return acceptableLibrarySpectrums;
//        }

//        private int GetFirstScanWithMassOverOrEqualInSpectralLibrary(double minimum)
//        {
//            int index = Array.BinarySearch(SpectralLibraryPrecursorMasses, minimum);
//            if (index < 0)
//            {
//                index = ~index;
//            }
//            // index of the first element that is larger than value
//            return index;
//        }


//    }
//}