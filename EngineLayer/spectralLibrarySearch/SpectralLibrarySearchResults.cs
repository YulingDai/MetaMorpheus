using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class SpectralLibrarySearchResults
    {

        public SpectralLibrarySearchResults(Spectrum theExperimentalSpectrum, int indexOfExperimentalSpectrum, List<SpectralLibrarayMatch> spectralLibrarayMatchs)
        {
            TheExperimentalSpectrum = theExperimentalSpectrum;
            IndexOfExperimentalSpectrum = indexOfExperimentalSpectrum;
            SpectralLibrarayMatchs = spectralLibrarayMatchs;
        }

        public Spectrum TheExperimentalSpectrum { get; set; }
        public int IndexOfExperimentalSpectrum { get; set; }
        public List<SpectralLibrarayMatch> SpectralLibrarayMatchs { get; set; }
}
}
