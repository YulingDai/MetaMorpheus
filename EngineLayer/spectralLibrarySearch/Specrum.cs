using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

namespace EngineLayer.spectralLibrarySearch
{
    public class Spectrum
    {
        public Spectrum()
        {
            //Name = name;
        }
        public Spectrum(String name, List< Peaks> peaks)
        {
            Name = name; ;
            PeaksInfo = peaks;
        }
        public string Name { get; set; }
        public double MW { get; set; }
        public double rententionTime { get; set; }
        public double precursorMz { get; set; }
        public int charge_state { get; set; }
        public List<Product> Peaks { get; set; }
        public Dictionary<Product, double> PeaksWithIntensity { get; set; }
        public double totalIonCurrent { get; set; }
        public double MonoisotopicMass { get; set; }
        public List<Peaks> PeaksInfo { get; set; }


        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Name + "\r\n");
            //spectrum.Append("precursor: " + precursorMz);
            spectrum.Append("Matched peaks number : " + this.PeaksInfo.Count + "\r\n");
            foreach (var eachPeak in this.PeaksInfo)
            {

                spectrum.Append(eachPeak.Mz + "\t" + eachPeak.Intensity + "\t" + "\"" + eachPeak.Product.ProductType.ToString() + eachPeak.Product.FragmentNumber.ToString() + "/" + 0 + "ppm" + "\"" + "\r\n");
            }

             return spectrum.ToString();

        }
    }
}
