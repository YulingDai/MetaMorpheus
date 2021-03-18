﻿using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer
{
    /// <summary>
    /// This class describes one element of a spectral library.
    /// </summary>
    public class LibrarySpectrum
    {

        public string Sequence { get; set; }
        public double RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }
        public string ModsString { get; set; }
        public string BaseSequenceWithoutMods { get; set; }

        public string Name
        {
            get { return Sequence + "/" + ChargeState; }
        }

        public LibrarySpectrum(string sequence, double precursorMz, int chargeState, List<MatchedFragmentIon> peaks, double rt, bool isDecoy = false)
        {
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            ChargeState = chargeState;
            IsDecoy = isDecoy;
            RetentionTime = rt;
        }

        //public override string ToString()
        //{
        //    StringBuilder spectrum = new StringBuilder();
        //    spectrum.Append("Name: " + Name);
        //    spectrum.Append("\nMW: " + PrecursorMz);
        //    spectrum.Append("\nComment: ");
        //    spectrum.Append("Parent=" + PrecursorMz);
        //    spectrum.Append(" RT=" + RetentionTime);
        //    spectrum.Append("\nNum peaks: " + MatchedFragmentIons.Count);

        //    double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();

        //    foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
        //    {
        //        double intensityFraction = matchedIon.Intensity / maxIntensity;

        //        spectrum.Append("\n" + matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
        //            matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
        //            matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
        //            matchedIon.Charge + "/" + 0 + "ppm" + "\"");
        //    }

        //    return spectrum.ToString();
        //}

        //test
        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Name);
            spectrum.Append("\nMW: " + PrecursorMz);
            spectrum.Append("\nComment: ");
            spectrum.Append("Parent=" + PrecursorMz);
            spectrum.Append(" RT=" + RetentionTime);
            spectrum.Append("\nNum peaks: " + MatchedFragmentIons.Count);

            double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();
            double sumIntensity = MatchedFragmentIons.Select(b => b.Intensity).Sum();

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                //double intensityFraction = matchedIon.Intensity / sumIntensity;
                double intensityFraction = matchedIon.Intensity ;

                spectrum.Append("\n" + matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + "/" + 0 + "ppm" + "\"");
            }

            return spectrum.ToString();
        }


    }
}