//using MassSpectrometry;
//using SharpLearning.Containers.Extensions;
//using System;
//using System.Collections.Generic;
//using System.Text;
//using ThermoFisher.CommonCore.Data.Interfaces;
//using System.Linq;
//using MzLibUtil;

//namespace EngineLayer.spectralLibrarySearch
//{
//    public class  ClassicSearchOfSpectralLibrary
//    {
//        private readonly double?[] SpectralLibraryPrecursorMasses;
//        public Spectrum[] sorted_spectralLibrary { get; set; }
//        public  SpectralLibrarySearchResults[] SpectralLibrarySearchResults { get; set; }
//        public double PrecursorTolerance { get; set; }
//        public double FragmentTolerance { get; set; }
//        public double SpectrumScoreCutOff { get; set; }
//        private readonly MassDiffAcceptor SearchMode;
//        public CommonParameters CommonParameters { get; set; }



//        public ClassicSearchOfSpectralLibrary(Spectrum[] spectralLibrary, Spectrum[] spectrums, double precursorTolerance, double fragmentTolerance, double spectrumScoreCutOff,
//            MassDiffAcceptor searchMode, CommonParameters commonParameters)
//        {
//            SearchMode = searchMode;
//            CommonParameters = commonParameters;
//            PrecursorTolerance = precursorTolerance;
//            FragmentTolerance = fragmentTolerance;
//            SpectrumScoreCutOff = spectrumScoreCutOff;
       
//            sorted_spectralLibrary = SortedByMzSpectralLibrary(spectralLibrary);
        
//            SpectralLibraryPrecursorMasses = sorted_spectralLibrary.Select(b => b.precursorMz).ToArray();
          
//            var spectrumsPrecursorMasses = spectrums.Select(b => b.precursorMz).ToArray();
            
//            SpectralLibrarySearchResults = new SpectralLibrarySearchResults[spectrums.Length];
          


//            for (int i = 0; i < spectrums.Length; i++)
//            {
//                var expecrimentalSpectrum = spectrums[i];
//                SpectralLibrarySearchResults[i] = new SpectralLibrarySearchResults(expecrimentalSpectrum, i, new List<SpectralLibrarayMatch>());
         
//                foreach (Spectrum matchedSpectrum in GetAcceptableLibrarySpectrums(Convert.ToDouble(expecrimentalSpectrum.precursorMz), SearchMode))
//                {
//                    List<PeaksInformationFromSpectrum> matchedPeaks = MatchSpectrumPeaks(expecrimentalSpectrum, matchedSpectrum);

//                    double thisScore = CalculateSpectrumScore(matchedSpectrum, matchedPeaks);
                  
//                    bool meetsScoreCutoff = thisScore >= 5;

//                    if (meetsScoreCutoff)
//                    {
       
//                        SpectralLibrarySearchResults[i].SpectralLibrarayMatchs.Add(new SpectralLibrarayMatch(matchedSpectrum, thisScore, matchedPeaks));
                
//                    }
    
//                }
        
//            }

//        }

//        //helper method for sorting
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
            
           
        
       

//        //public Spectrum[] MatchedSpectrums(Spectrum toSearch, Spectrum[] sorted_SpectralLibrary)
//        //{
//        //    var spectrums = new List<Spectrum>();
//        //    foreach(Spectrum x in sorted_SpectralLibrary)
//        //    {
//        //        if (Math.Abs(Convert.ToDouble(toSearch.precursorMz) - Convert.ToDouble(x.precursorMz)) <= PrecursorTolerance)
//        //        {
//        //            spectrums.Add(x);
//        //        }
//        //    }
//        //    return spectrums.ToArray();
//        //}

//        public static double CalculateSpectrumScore(Spectrum thisSpectrum, List<PeaksInformationFromSpectrum> matchedPeaks)
//        {
//            double score = 0;
//            for (int i = 0; i < matchedPeaks.Count; i++)
//            {
//                score += 1 + matchedPeaks[i].Intensity;
//            }
//            return score;
//        }

//        //public List<Spectrum> GetAcceptableLibrarySpectrums(double experimentalPreMz, double massError)
//        //{
//        //    var acceptableLibrarySpectrums = new List<Spectrum>();
//        //    var allowedPreMz = new DoubleRange(experimentalPreMz - massError, experimentalPreMz + massError);
//        //    int index = GetFirstScanWithMassOverOrEqualInSpectralLibrary(allowedPreMz.Minimum);
     
//        //    if (index< SpectralLibraryPrecursorMasses.Length)
//        //    {
//        //        var precursorMass = SpectralLibraryPrecursorMasses[index];
//        //        while (precursorMass <= allowedPreMz.Maximum)
//        //        {
//        //            Spectrum matchedSpectrum = sorted_spectralLibrary[index];
//        //            acceptableLibrarySpectrums.Add(matchedSpectrum);
//        //            index++;
//        //            if (index == sorted_spectralLibrary.Length)
//        //            {
//        //                break;
//        //            }

//        //            precursorMass = SpectralLibraryPrecursorMasses[index];
//        //        }
//        //    }
//        //    return acceptableLibrarySpectrums;
//        //}

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
//                        Spectrum matchedSpectrum = sorted_spectralLibrary[index];
//                        acceptableLibrarySpectrums.Add(matchedSpectrum);
//                        index++;
//                        if (index == sorted_spectralLibrary.Length)
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
