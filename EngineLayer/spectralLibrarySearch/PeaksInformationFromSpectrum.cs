using System;
using System.Collections.Generic;
using System.Text;
using Proteomics.Fragmentation;

namespace EngineLayer.spectralLibrarySearch
{
    public class Peaks
    {
        public double Mz { get; set; }
        public double Intensity { get; set; }
        public string SpectrumPeakProductType { get; set; }
        public int FragmentNumber { get; set; }
        public double massErrorPpm { get; set; }
        public readonly int AminoAcidPosition;
        public Product Product { get; set; }

        public Peaks(double mz, double intensity, Product product)
        {
            Mz = mz;
            Intensity = intensity;
            Product = product;
        }

        //public override bool Equals(object obj);
        //public override int GetHashCode();
        //public override string ToString();

    }
}