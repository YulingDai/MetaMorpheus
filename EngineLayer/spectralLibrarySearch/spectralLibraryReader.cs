using Easy.Common.Extensions;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.spectralLibrarySearch
{
    public class spectralLibraryReader
    {
        public Spectrum[] spectrums { get; set; }

        public spectralLibraryReader(string filePath)
        {
            var spectrumsList = new List<Spectrum>();
            string[] lines;

            
            try
            {
                lines = File.ReadAllLines(filePath);
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not read file: " + e.Message);
            }

            //find the lines which contain "name"
            var nameLine = new List<int>();
            for (int i = 0; i < lines.Length; i++)
            {
                if (lines[i].Contains("name", StringComparison.OrdinalIgnoreCase))
                {
                    nameLine.Add(i);
                }
            }
            nameLine.Add(lines.Length);// for the convenience to separate the file to different parts
            //nameLine.ForEach(i => Console.Write("{0}\t", i));

            //load separate spectrums 
            for (int i = 0; i < nameLine.Count-1; i++)
            {
                //load each spectrum
                var singleSpectrum = new Spectrum();
                for (int j = nameLine[i]; j < nameLine[i + 1]; j++)
                {

                    //get name of each spectrum
                    if (lines[j].Contains("name", StringComparison.OrdinalIgnoreCase))
                    {
                        try 
                        {
                            string[] name = lines[j].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                            singleSpectrum.Name = name[1];
                        }
                        catch(Exception e)
                        {
                            singleSpectrum.Name = null;
                            throw new MetaMorpheusException("Could not find the name : " + e.Message);
                        }
                    }

                    //get MW of each spectrum
                    if (lines[j].Contains("MW", StringComparison.OrdinalIgnoreCase) || lines[j].Contains("Monoisotopic Mass", StringComparison.OrdinalIgnoreCase))
                    {
                        try
                        {
                            string[] mw = lines[j].Split(":").Select(b => b.Trim()).ToArray();
                            singleSpectrum.MW = double.Parse(mw[1]);
                        }
                        catch(Exception e)
                        {
                            singleSpectrum.MW = null;
                            throw new MetaMorpheusException("Could not find the MW : " + e.Message);
                        }
                    }

                    // get information from comment
                    if (lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                    {
                     
                        string[] comment = lines[j].Split(" ").Select(b => b.Trim()).ToArray();
                        for (int l = 0; l< comment.Length; l++)
                        {
                            try
                            {
                                if (comment[l].Contains("charge", StringComparison.OrdinalIgnoreCase))
                                {
                                    string[] charge_state = comment[l].Split(":").Select(b => b.Trim()).ToArray();
                                    singleSpectrum.charge_state = int.Parse(charge_state[1]);
                                }
                            }
                            catch(Exception e)
                            {
                                singleSpectrum.charge_state = null;
                                throw new MetaMorpheusException("Could not find the charge state : " + l + e.Message);
                            }

                            try
                            {
                                if (comment[l].Contains("parent", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("precursor", StringComparison.OrdinalIgnoreCase))
                                {
                                    string[] precursorMz = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                    singleSpectrum.precursorMz = double.Parse(precursorMz[1]);
                                }
                            }
                            catch (Exception e)
                            {
                                singleSpectrum.precursorMz = null;
                                throw new MetaMorpheusException("Could not find the mz of precursor : " + e.Message);
                            }

                            try
                            {
                                if (comment[l].Contains("RT", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("retention time", StringComparison.OrdinalIgnoreCase))
                                {
                                    string[] rententionTime = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                    singleSpectrum.rententionTime = double.Parse(rententionTime[1]);
                                }
                            }
                            catch (Exception e)
                            {
                                singleSpectrum.rententionTime = null;
                                throw new MetaMorpheusException("Could not find the mz of rentention time : " + e.Message);
                            }
                        }
                    }

                    if (lines[j].Contains("peaks", StringComparison.OrdinalIgnoreCase))
                    {
                        string[] numPeaks = lines[j].Split(":").Select(b => b.Trim()).ToArray();
                        var numberOfPeaks = int.Parse(numPeaks[1]);
                        var peaksList = new List<PeaksInformationFromSpectrum>();
                        for (int k = j+1; k < j+1+numberOfPeaks; k++)
                        {
                            string[] eachPeak = lines[k].Split("\t").Select(b => b.Trim()).ToArray();
                            var b = new PeaksInformationFromSpectrum(double.Parse(eachPeak[0]), double.Parse(eachPeak[1]));
                 
                            string[] by = eachPeak[2].Split(new char[] { '/', '\"', 'p'}, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                            b.massErrorPpm = double.Parse(by[1].Trim());
                            peaksList.Add(b);
                        }
                        singleSpectrum.Peaks = peaksList.ToArray();
                    }
                }
                if(singleSpectrum.Name!=null && singleSpectrum.precursorMz!=null&&singleSpectrum.Peaks.Length!=0)
                {
                    spectrumsList.Add(singleSpectrum);
                }
            }
            if(spectrumsList.Count==0)
            {
                this.spectrums = null;
                Console.WriteLine("the library doesn't contain any spectrums!");
            }
            else
            {
                this.spectrums = spectrumsList.ToArray();
            }
        }
    }
}
