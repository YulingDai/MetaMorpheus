﻿using Chemistry;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {
        public Ms2ScanWithSpecificMass(MsDataScan mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath, CommonParameters commonParam, IsotopicEnvelope[] neutralExperimentalFragments = null)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            FullFilePath = fullFilePath;
            ChildScans = new List<Ms2ScanWithSpecificMass>();
            NativeId = mzLibScan.NativeId;

            TheScan = mzLibScan;

            if (commonParam.DissociationType != DissociationType.LowCID)
            {
                ExperimentalFragments = neutralExperimentalFragments ?? GetNeutralExperimentalFragments(mzLibScan, commonParam);
            }

            if (ExperimentalFragments != null && ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
            }
            else
            {
                DeconvolutedMonoisotopicMasses = new double[0];
            }
        }

        public MsDataScan TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        public List<Ms2ScanWithSpecificMass> ChildScans { get; set; } // MS2/MS3 scans that are children of this MS2 scan
        private double[] DeconvolutedMonoisotopicMasses;
        public string NativeId { get; } 

        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;

        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;

        public double RetentionTime => TheScan.RetentionTime;

        public int NumPeaks => TheScan.MassSpectrum.Size;

        public double TotalIonCurrent => TheScan.TotalIonCurrent;

        public static IsotopicEnvelope[] GetNeutralExperimentalFragments(MsDataScan scan, CommonParameters commonParam)
        {
            int minZ = 1;
            int maxZ = 10;

            var neutralExperimentalFragmentMasses = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range,
                minZ, maxZ, commonParam.DeconvolutionMassTolerance.Value, commonParam.DeconvolutionIntensityRatio).ToList();

            if (commonParam.AssumeOrphanPeaksAreZ1Fragments)
            {
                HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                    .SelectMany(p => p.Peaks.Select(v => ClassExtensions.RoundedDouble(v.mz).Value)));

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    double mz = scan.MassSpectrum.XArray[i];
                    double intensity = scan.MassSpectrum.YArray[i];

                    if (!alreadyClaimedMzs.Contains(ClassExtensions.RoundedDouble(mz).Value))
                    {
                        neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (mz, intensity) },
                            mz.ToMass(1), 1, intensity, 0, 0));
                    }
                }
            }

            return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        public IsotopicEnvelope GetClosestExperimentalIsotopicEnvelope(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass)];
        }
        public IsotopicEnvelope GetClosestExperimentalIsotopicEnvelopeWithChargeState1(double theoreticalNeutralMass)
        {
            var Charge1_ExperimentalFragments = ExperimentalFragments.ToList().Where(p => p.Charge == 1).ToArray();
            double[] Charge1_DeconvolutedMonoisotopicMasses = Charge1_ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
            if (Charge1_DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return Charge1_ExperimentalFragments[GetClosestFragmentMass_charge1(theoreticalNeutralMass, Charge1_DeconvolutedMonoisotopicMasses)];
        }
        public int GetClosestFragmentMass(double mass)
        {
            int index = Array.BinarySearch(DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index == DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }
        public int GetClosestFragmentMass_charge1(double mass, double[] Charge1_DeconvolutedMonoisotopicMasses)
        {
            int index = Array.BinarySearch(Charge1_DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index == Charge1_DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - Charge1_DeconvolutedMonoisotopicMasses[index - 1] > Charge1_DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }
        public double? GetClosestExperimentalFragmentMz(double theoreticalMz, out double? intensity)
        {
            if (TheScan.MassSpectrum.XArray.Length == 0)
            {
                intensity = null;
                return null;
            }
            intensity = TheScan.MassSpectrum.YArray[GetClosestFragmentMzIndex(theoreticalMz).Value];
            return TheScan.MassSpectrum.XArray[GetClosestFragmentMzIndex(theoreticalMz).Value];
        }

        private int? GetClosestFragmentMzIndex(double mz)
        {
            if (TheScan.MassSpectrum.XArray.Length == 0)
            {
                return null;
            }
            int index = Array.BinarySearch(TheScan.MassSpectrum.XArray, mz);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index >= TheScan.MassSpectrum.XArray.Length)
            {
                return index - 1;
            }
            if (index == 0)
            {
                return index;
            }

            if (mz - TheScan.MassSpectrum.XArray[index - 1] > TheScan.MassSpectrum.XArray[index] - mz)
            {
                return index;
            }
            return index - 1;

        }
    }
}