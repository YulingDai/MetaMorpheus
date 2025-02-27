﻿using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        public SearchTask() : base(MyTask.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();
        }

        public SearchParameters SearchParameters { get; set; }

        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                    {
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    }
                    else
                    {
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);
                    }

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                case MassDiffAcceptorType.PlusOrMinusThreeMM:
                    return new DotMassDiffAcceptor(
                        "PlusOrMinus3Da",
                        new List<double>
                        {
                            -3 * Chemistry.Constants.C13MinusC12,
                            -2 * Chemistry.Constants.C13MinusC12,
                            -1 * Chemistry.Constants.C13MinusC12,
                            0,
                            1 * Chemistry.Constants.C13MinusC12,
                            2 * Chemistry.Constants.C13MinusC12,
                            3 * Chemistry.Constants.C13MinusC12
                        },
                        precursorMassTolerance);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            if (SearchParameters.DoQuantification)
            {
                // disable quantification if a .mgf is being used
                if (currentRawFileList.Any(x => Path.GetExtension(x).Equals(".mgf", StringComparison.OrdinalIgnoreCase)))
                {
                    SearchParameters.DoQuantification = false;
                }
                //if we're doing SILAC, assign and add the silac labels to the residue dictionary
                else if (SearchParameters.SilacLabels != null || SearchParameters.StartTurnoverLabel != null || SearchParameters.EndTurnoverLabel != null)
                {
                    char heavyLabel = 'a'; //char to assign
                    //add the Turnoverlabels to the silacLabels list. They weren't there before just to prevent duplication in the tomls
                    if (SearchParameters.StartTurnoverLabel != null || SearchParameters.EndTurnoverLabel != null)
                    {
                        //original silacLabels object is null, so we need to initialize it
                        SearchParameters.SilacLabels = new List<SilacLabel>();
                        if (SearchParameters.StartTurnoverLabel != null)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.StartTurnoverLabel, heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            SearchParameters.StartTurnoverLabel = updatedLabel.updatedLabel;
                            SearchParameters.SilacLabels.Add(SearchParameters.StartTurnoverLabel);
                        }
                        if (SearchParameters.EndTurnoverLabel != null)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.EndTurnoverLabel, heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            SearchParameters.EndTurnoverLabel = updatedLabel.updatedLabel;
                            SearchParameters.SilacLabels.Add(SearchParameters.EndTurnoverLabel);
                        }
                    }
                    else
                    {
                        //change the silac residues to lower case amino acids (currently null)
                        List<SilacLabel> updatedLabels = new List<SilacLabel>();
                        for (int i = 0; i < SearchParameters.SilacLabels.Count; i++)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.SilacLabels[i], heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            updatedLabels.Add(updatedLabel.updatedLabel);
                        }
                        SearchParameters.SilacLabels = updatedLabels;
                    }
                }
            }
            //if no quant, remove any silac labels that may have been added, because they screw up downstream analysis
            if (!SearchParameters.DoQuantification) //using "if" instead of "else", because DoQuantification can change if it's an mgf
            {
                SearchParameters.SilacLabels = null;
            }

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);
            SanitizeProteinDatabase(proteinList, SearchParameters.TCAmbiguity);

            // load spectral libraries
            var spectralLibrary = LoadSpectralLibraries(taskId, dbFilenameList);

            // write prose settings
            ProseCreatedWhileRunning.Append("The following search settings were used: ");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("report PSM ambiguity = " + CommonParameters.ReportAllAmbiguity + ". ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy)
                + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the search task
            MyTaskResults = new MyTaskResults(this);
            List<PeptideSpectralMatch> allPsms = new List<PeptideSpectralMatch>();

            //generate an array to store category specific fdr values (for speedy semi/nonspecific searches)
            int numFdrCategories = (int)(Enum.GetValues(typeof(FdrCategory)).Cast<FdrCategory>().Last() + 1); //+1 because it starts at zero
            List<PeptideSpectralMatch>[] allCategorySpecificPsms = new List<PeptideSpectralMatch>[numFdrCategories];
            for (int i = 0; i < numFdrCategories; i++)
            {
                allCategorySpecificPsms[i] = new List<PeptideSpectralMatch>();
            }

            FlashLfqResults flashLfqResults = null;

            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.PrecursorMass).ToArray();
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                PeptideSpectralMatch[] fileSpecificPsms = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // modern search
                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / combinedParams.TotalPartitions,
                            ((currentPartition + 1) * proteinList.Count / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count / combinedParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                            SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, combinedParams, FileSpecificParameters,
                            SearchParameters.MaxFragmentSize, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        List<int>[] precursorIndex = null;

                        lock (indexLock)
                        {
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                        }

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition,
                            combinedParams, this.FileSpecificParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }
                }
                // nonspecific search
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    PeptideSpectralMatch[][] fileSpecificPsmsSeparatedByFdrCategory = new PeptideSpectralMatch[numFdrCategories][]; //generate an array of all possible locals
                    for (int i = 0; i < numFdrCategories; i++) //only add if we're using for FDR, else ignore it as null.
                    {
                        fileSpecificPsmsSeparatedByFdrCategory[i] = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                    }

                    //create params for N, C, or both if semi
                    List<CommonParameters> paramsToUse = new List<CommonParameters> { combinedParams };
                    if (combinedParams.DigestionParams.SearchModeType == CleavageSpecificity.Semi) //if semi, we need to do both N and C to hit everything
                    {
                        paramsToUse.Clear();
                        List<FragmentationTerminus> terminiToUse = new List<FragmentationTerminus> { FragmentationTerminus.N, FragmentationTerminus.C };
                        foreach (FragmentationTerminus terminus in terminiToUse) //set both termini
                        {
                            paramsToUse.Add(combinedParams.CloneWithNewTerminus(terminus));
                        }
                    }

                    //Compress array of deconvoluted ms2 scans to avoid searching the same ms2 multiple times while still identifying coisolated peptides
                    List<int>[] coisolationIndex = new List<int>[] { new List<int>() };
                    if (arrayOfMs2ScansSortedByMass.Length != 0)
                    {
                        int maxScanNumber = arrayOfMs2ScansSortedByMass.Max(x => x.OneBasedScanNumber);
                        coisolationIndex = new List<int>[maxScanNumber + 1];
                        for (int i = 0; i < arrayOfMs2ScansSortedByMass.Length; i++)
                        {
                            int scanNumber = arrayOfMs2ScansSortedByMass[i].OneBasedScanNumber;
                            if (coisolationIndex[scanNumber] == null)
                            {
                                coisolationIndex[scanNumber] = new List<int> { i };
                            }
                            else
                            {
                                coisolationIndex[scanNumber].Add(i);
                            }
                        }
                        coisolationIndex = coisolationIndex.Where(x => x != null).ToArray();
                    }

                    //foreach terminus we're going to look at
                    foreach (CommonParameters paramToUse in paramsToUse)
                    {
                        //foreach database partition
                        for (int currentPartition = 0; currentPartition < paramToUse.TotalPartitions; currentPartition++)
                        {
                            List<PeptideWithSetModifications> peptideIndex = null;

                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / paramToUse.TotalPartitions,
                                ((currentPartition + 1) * proteinList.Count / paramToUse.TotalPartitions) - (currentPartition * proteinList.Count / paramToUse.TotalPartitions));

                            List<int>[] fragmentIndex = null;
                            List<int>[] precursorIndex = null;

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                                SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, paramToUse, FileSpecificParameters,
                                SearchParameters.MaxFragmentSize, true, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                            lock (indexLock)
                            {
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                            }

                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsmsSeparatedByFdrCategory, arrayOfMs2ScansSortedByMass, coisolationIndex, peptideIndex, fragmentIndex,
                                precursorIndex, currentPartition, paramToUse, this.FileSpecificParameters, variableModifications, massDiffAcceptor,
                                SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + paramToUse.TotalPartitions + "!", thisId));
                            if (GlobalVariables.StopLoops) { break; }
                        }
                    }
                    lock (psmLock)
                    {
                        for (int i = 0; i < allCategorySpecificPsms.Length; i++)
                        {
                            if (allCategorySpecificPsms[i] != null)
                            {
                                allCategorySpecificPsms[i].AddRange(fileSpecificPsmsSeparatedByFdrCategory[i]);
                            }
                        }
                    }
                }
                // classic search
                else
                {
                    Status("Starting search...", thisId);
                    var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                       SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, proteinList, massDiffAcceptor, combinedParams, this.FileSpecificParameters, spectralLibrary, thisId,SearchParameters.WriteSpectralLibrary);
                    newClassicSearchEngine.Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                //look for internal fragments
                if (SearchParameters.MinAllowedInternalFragmentLength != 0)
                {
                    MatchInternalFragmentIons(fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams, SearchParameters.MinAllowedInternalFragmentLength);
                }

                // calculate/set spectral angles if there is a spectral library being used
                if (spectralLibrary != null)
                {
                    Status("Calculating spectral library similarity...", thisId);
                }
                SpectralLibrarySearchFunction.CalculateSpectralAngles(spectralLibrary, fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams);

                lock (psmLock)
                {
                    allPsms.AddRange(fileSpecificPsms);
                }

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            if (spectralLibrary != null)
            {
                spectralLibrary.CloseConnections();
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            int numNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            //resolve category specific fdrs (for speedy semi and nonspecific
            if (SearchParameters.SearchType == SearchType.NonSpecific)
            {
                allPsms = NonSpecificEnzymeSearchEngine.ResolveFdrCategorySpecificPsms(allCategorySpecificPsms, numNotches, taskId, CommonParameters, FileSpecificParameters);
            }

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters
            {
                SearchTaskResults = MyTaskResults,
                SearchTaskId = taskId,
                SearchParameters = SearchParameters,
                ProteinList = proteinList,
                AllPsms = allPsms,
                VariableModifications = variableModifications,
                FixedModifications = fixedModifications,
                ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams)),
                CurrentRawFileList = currentRawFileList,
                MyFileManager = myFileManager,
                NumNotches = numNotches,
                OutputFolder = OutputFolder,
                IndividualResultsOutputFolder = Path.Combine(OutputFolder, "Individual File Results"),
                FlashLfqResults = flashLfqResults,
                FileSettingsList = fileSettingsList,
                NumMs2SpectraPerFile = numMs2SpectraPerFile,
                DatabaseFilenameList = dbFilenameList
            };
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask
            {
                Parameters = parameters,
                FileSpecificParameters = this.FileSpecificParameters,
                CommonParameters = CommonParameters
            };
            return postProcessing.Run();
        }

        private int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact: return 1;
                case MassDiffAcceptorType.OneMM: return 2;
                case MassDiffAcceptorType.TwoMM: return 3;
                case MassDiffAcceptorType.ThreeMM: return 4;
                case MassDiffAcceptorType.ModOpen: return 1;
                case MassDiffAcceptorType.Open: return 1;
                case MassDiffAcceptorType.PlusOrMinusThreeMM: return 7;
                case MassDiffAcceptorType.Custom: return ParseSearchMode(customMdac).NumNotches;

                default: throw new MetaMorpheusException("Unknown mass difference acceptor type");
            }
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor massDiffAcceptor = null;

            try
            {
                var split = text.Split(' ');

                switch (split[1])
                {
                    case "dot":
                        double[] massShifts = split[4].Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToArray();
                        string newString = split[2].Replace("�", "");
                        double toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                        if (split[3].ToUpperInvariant().Equals("PPM"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                        }
                        else if (split[3].ToUpperInvariant().Equals("DA"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                        }

                        break;

                    case "interval":
                        IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(';'), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[0],
                            CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[1], CultureInfo.InvariantCulture)));
                        massDiffAcceptor = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                        break;

                    case "OpenSearch":
                        massDiffAcceptor = new OpenSearchMode();
                        break;

                    case "daltonsAroundZero":
                        massDiffAcceptor = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    case "ppmAroundZero":
                        massDiffAcceptor = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    default:
                        throw new MetaMorpheusException("Unrecognized search mode type: " + split[1]);
                }
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not parse search mode string: " + e.Message);
            }

            return massDiffAcceptor;
        }

        public static void MatchInternalFragmentIons(PeptideSpectralMatch[] fileSpecificPsms, Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass, CommonParameters combinedParams, int minInternalFragmentLength)
        {
            //for each PSM with an ID
            for (int index = 0; index < fileSpecificPsms.Length; index++)
            {
                PeptideSpectralMatch psm = fileSpecificPsms[index];
                if (psm != null && psm.BestMatchingPeptides.Count() > 0)
                {
                    //Get the scan
                    Ms2ScanWithSpecificMass scanForThisPsm = arrayOfMs2ScansSortedByMass[index];
                    DissociationType dissociationType = combinedParams.DissociationType == DissociationType.Autodetect ?
                    scanForThisPsm.TheScan.DissociationType.Value : combinedParams.DissociationType;

                    //Get the theoretical peptides
                    List<PeptideWithSetModifications> ambiguousPeptides = new List<PeptideWithSetModifications>();
                    List<int> notches = new List<int>();
                    foreach (var (Notch, Peptide) in psm.BestMatchingPeptides)
                    {
                        ambiguousPeptides.Add(Peptide);
                        notches.Add(Notch);
                    }

                    //get matched ions for each peptide
                    List<List<MatchedFragmentIon>> matchedIonsForAllAmbiguousPeptides = new List<List<MatchedFragmentIon>>();
                    List<Product> internalFragments = new List<Product>();
                    foreach (PeptideWithSetModifications peptide in ambiguousPeptides)
                    {
                        internalFragments.Clear();
                        peptide.FragmentInternally(combinedParams.DissociationType, minInternalFragmentLength, internalFragments);
                        //TODO: currently, internal and terminal ions can match to the same observed peaks (much like how b- and y-ions can match to the same peaks). Investigate if we should change that...                        
                        matchedIonsForAllAmbiguousPeptides.Add(MetaMorpheusEngine.MatchFragmentIons(scanForThisPsm, internalFragments, combinedParams));
                    }

                    //Find the max number of matched ions
                    int maxNumMatchedIons = matchedIonsForAllAmbiguousPeptides.Max(x => x.Count);

                    //remove peptides if they have fewer than max-1 matched ions, thus requiring at least two internal ions to disambiguate an ID
                    //if not removed, then add the matched internal ions
                    HashSet<PeptideWithSetModifications> PeptidesToMatchingInternalFragments = new HashSet<PeptideWithSetModifications>();
                    for (int peptideIndex = 0; peptideIndex < ambiguousPeptides.Count; peptideIndex++)
                    {
                        //if we should remove the theoretical, remove it
                        if (matchedIonsForAllAmbiguousPeptides[peptideIndex].Count + 1 < maxNumMatchedIons)
                        {
                            psm.RemoveThisAmbiguousPeptide(notches[peptideIndex], ambiguousPeptides[peptideIndex]);
                        }
                        // otherwise add the matched internal ions to the total ions
                        else
                        {
                            PeptideWithSetModifications currentPwsm = ambiguousPeptides[peptideIndex];
                            //check that we haven't already added the matched ions for this peptide
                            if (!PeptidesToMatchingInternalFragments.Contains(currentPwsm))
                            {
                                PeptidesToMatchingInternalFragments.Add(currentPwsm); //record that we've seen this peptide
                                psm.PeptidesToMatchingFragments[currentPwsm].AddRange(matchedIonsForAllAmbiguousPeptides[peptideIndex]); //add the matched ions
                            }
                        }
                    }
                }
            }
        }
    }
}