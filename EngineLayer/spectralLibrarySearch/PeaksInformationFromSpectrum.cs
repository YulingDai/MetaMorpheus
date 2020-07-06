using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class PeaksInformationFromSpectrum
    {
        public double Mz { get; set; }
        public double Intensity { get; set; }
        public string  SpectrumPeakProductType { get; set; }
        public int fragmentNumber { get; set; }
        public double massErrorPpm { get; set; }

        public PeaksInformationFromSpectrum(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

    }
}
