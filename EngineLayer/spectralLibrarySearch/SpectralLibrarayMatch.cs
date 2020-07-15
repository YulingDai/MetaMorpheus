using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectralLibrarayMatch
    {

        //public Spectrum TheExperimentalSpectrum { get; set; }
        //public int IndexOfExperimentalSpectrum { get; set; }
        public Spectrum MatchedSpectrumFromLibrary { get; set; }
        public double MatchScore { get; set; }
        public List<MatchedFragmentIon> MatchedPeaks { get; set; }
      
        public SpectralLibrarayMatch( Spectrum matchedSpectrumFromLibrary, double matchScore, List<MatchedFragmentIon> matchedPeaks)
        {
            //TheExperimentalSpectrum = theExperimentalSpectrum;
            //IndexOfExperimentalSpectrum = indexOfExperimentalSpectrum;
            MatchedSpectrumFromLibrary = matchedSpectrumFromLibrary;
            MatchScore = matchScore;
            MatchedPeaks = matchedPeaks;

        }
    }
}
