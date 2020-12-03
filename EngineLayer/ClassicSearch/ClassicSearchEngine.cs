﻿using Chemistry;
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
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly List<Protein> Proteins;
        private readonly List<Modification> FixedModifications;
        private readonly List<Modification> VariableModifications;
        private readonly List<SilacLabel> SilacLabels;
        private readonly (SilacLabel StartLabel, SilacLabel EndLabel)? TurnoverLabels;
        private readonly PeptideSpectralMatch[] PeptideSpectralMatches;
        private readonly Ms2ScanWithSpecificMass[] ArrayOfSortedMS2Scans;
        private readonly double[] MyScanPrecursorMasses;

        public ClassicSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans,
            List<Modification> variableModifications, List<Modification> fixedModifications, List<SilacLabel> silacLabels, SilacLabel startLabel, SilacLabel endLabel,
            List<Protein> proteinList, MassDiffAcceptor searchMode, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            SpectralLibrary spectralLibrary, List<string> nestedIds)
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
            SpectralLibrary = spectralLibrary;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
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
                                //if (SpectralLibrary != null && !SpectralLibrary.ContainsSpectrum(peptide.FullSequence, scan.TheScan.PrecursorCharge))
                                //{
                                    //continue;
                                //}

                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideTheorProducts, CommonParameters);

                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);
                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

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
                                        }
                                    }
                                }
                            }

                            List<Product> peptideDecoyProducts = GetDecoyProductsNB(peptideTheorProducts);
                            Protein decoyProtein = new Protein(peptide.Protein.BaseSequence, "DECOY_" + peptide.Protein.Accession, null, null, null, null, null, null, false, false, null, null, null, null, null, null, "");
                            PeptideWithSetModifications decoyPeptide = new PeptideWithSetModifications(decoyProtein, new DigestionParams(), peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, CleavageSpecificity.Full, peptide.PeptideDescription, peptide.MissedCleavages, new Dictionary<int, Modification>(), peptide.NumFixedMods, null);

                            foreach (ScanWithIndexAndNotchInfo scan in GetAcceptableScans(peptide.MonoisotopicMass, SearchMode))
                            {
                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan, peptideDecoyProducts, CommonParameters);

                                double thisScore = CalculatePeptideScore(scan.TheScan.TheScan, matchedIons);

                                bool meetsScoreCutoff = thisScore >= CommonParameters.ScoreCutoff;

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
                                                PeptideSpectralMatches[scan.ScanIndex] = new PeptideSpectralMatch(decoyPeptide, scan.Notch, thisScore, scan.ScanIndex, scan.TheScan, CommonParameters.DigestionParams, matchedIons, 0);
                                            }
                                            else
                                            {
                                                PeptideSpectralMatches[scan.ScanIndex].AddOrReplace(decoyPeptide, thisScore, scan.Notch, CommonParameters.ReportAllAmbiguity, matchedIons, 0);
                                            }
                                        }
                                    }
                                }
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

            foreach (PeptideSpectralMatch psm in PeptideSpectralMatches.Where(p => p != null))
            {
                psm.ResolveAllAmbiguities();
            }

            return new MetaMorpheusEngineResults(this);
        }

         private List<Product> GetDecoyProductsNB(List<Product> products)
        {
            List<Product> N_terminalTheoreticalProducts = products.Where(p => p.Terminus == FragmentationTerminus.N).ToList();

            List<double> masses = N_terminalTheoreticalProducts.Select(m => m.NeutralMass).ToList();
            masses = masses.OrderBy(m => m).ToList();
            List<double> massDifferences = new List<double>();
            for (int i = 0; i < masses.Count-1; i++)
            {
                massDifferences.Add(masses[i + 1] - masses[i]);
            }

            massDifferences = ShuffleList(massDifferences);

            List<double> newMasses = new List<double>();
            newMasses.Add(masses.Max());
            double firstMass = masses.Max();
            foreach (double massDifference in massDifferences)
            {
                firstMass -= massDifference;
                newMasses.Add(firstMass);
            }

            List<Product> newProducts = new List<Product>();
            for (int i = 0; i < N_terminalTheoreticalProducts.Count; i++)
            {
                Product original = N_terminalTheoreticalProducts[i];
                Product newProduct = new Product(original.ProductType, original.Terminus, newMasses[i], original.FragmentNumber, original.AminoAcidPosition, original.NeutralLoss);
                newProducts.Add(newProduct);
            }

            newProducts.AddRange(GetDecoyProductsCY(products.Where(p => p.Terminus == FragmentationTerminus.C).ToList(), massDifferences));

            return newProducts;
        }

        private List<Product> GetDecoyProductsCY(List<Product> peptideTheorProducts, List<double> massDifferences)
        {
            List<double> masses = peptideTheorProducts.Select(m => m.NeutralMass).ToList();
            masses = masses.OrderBy(m => m).ToList();

            List<double> newMasses = new List<double>();

            double firstMass = masses.Min() + ChemicalFormula.ParseFormula("H-2O-1").MonoisotopicMass;
            newMasses.Add(firstMass);

            foreach (double massDifference in massDifferences)
            {
                firstMass += (massDifference);
                newMasses.Add(firstMass);
            }

            List<Product> newProducts = new List<Product>();
            for (int i = 0; i < peptideTheorProducts.Count; i++)
            {
                Product original = peptideTheorProducts[i];
                Product newProduct = new Product(original.ProductType, original.Terminus, newMasses[i], original.FragmentNumber, original.AminoAcidPosition, original.NeutralLoss);
                newProducts.Add(newProduct);
            }
            return newProducts;
        }

        private List<E> ShuffleList<E>(List<E> inputList)
        {
            List<E> randomList = new List<E>();

            Random r = new Random();
            int randomIndex = 0;
            while (inputList.Count > 0)
            {
                randomIndex = r.Next(0, inputList.Count); //Choose a random object in the list
                randomList.Add(inputList[randomIndex]); //add it to the new, random list
                inputList.RemoveAt(randomIndex); //remove to avoid duplicates
            }

            return randomList; //return the new random list
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