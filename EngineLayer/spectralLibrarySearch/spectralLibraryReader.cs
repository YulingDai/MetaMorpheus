using Easy.Common.Extensions;
using Proteomics.Fragmentation;
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
        public Dictionary<String, Spectrum> SpectralLibraryDictionary { get; set; }
        //int i = 1;

        public spectralLibraryReader(string filePath)
        {
            SpectralLibraryDictionary = new Dictionary<string, Spectrum>();
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
            Console.WriteLine(nameLine.Count);

            //load separate spectrums 
            for (int i = 0; i < nameLine.Count - 1; i++)
                {
                    //load each spectrum
                    var singleSpectrum = new Spectrum();
                    for (int j = nameLine[i]; j < nameLine[i + 1]-1; j++)
                    {
                        //get name of each spectrum
                        if (lines[j].Contains("name", StringComparison.OrdinalIgnoreCase))
                        {
                            try
                            {
                                string[] name = lines[j].Split(new char[] { ':', '=' }, 2).Select(b => b.Trim()).ToArray();
                                singleSpectrum.Name = name[1];
                            }
                            catch (Exception e)
                            {
                                //singleSpectrum.Name = null;
                                throw new MetaMorpheusException("Could not find the name : " + e.Message);
                            }
                        }

                        //get MW of each spectrum
                        else if ((lines[j].Contains("MW", StringComparison.OrdinalIgnoreCase) || lines[j].Contains("Monoisotopic Mass", StringComparison.OrdinalIgnoreCase))&& !lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                        {
                        string[] mw = lines[j].Split(":", 2).Select(b => b.Trim()).ToArray();
                            try
                            {
                                var x = mw[0] + mw[1];
                                singleSpectrum.MW = double.Parse(mw[1]);
                            }
                            catch (Exception e)
                            {
                            //singleSpectrum.MW = 0;
                            //Console.WriteLine(lines[j]+ "  " + mw[0] + "  " + mw[1]);
                                throw new MetaMorpheusException("Could not find the MW : " + e.Message);
                            }
                        }

                        // get information from comment
                        if (lines[j].Contains("comment", StringComparison.OrdinalIgnoreCase))
                        {

                            string[] comment = lines[j].Split(" ").Select(b => b.Trim()).ToArray();
                            for (int l = 0; l < comment.Length; l++)
                            {
                                if (comment[l].Contains("charge", StringComparison.OrdinalIgnoreCase))
                                {
                                    try
                                    {
                                        string[] charge_state = comment[l].Split(":").Select(b => b.Trim()).ToArray();
                                        singleSpectrum.charge_state = int.Parse(charge_state[1]);
                                    }
                                    catch (Exception e)
                                    {
                                        //singleSpectrum.charge_state = 0;
                                        throw new MetaMorpheusException("Could not find the charge state : " + e.Message);
                                    }
                                }

                                if (comment[l].Contains("parent", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("precursor", StringComparison.OrdinalIgnoreCase))
                                {
                                    try
                                    {
                                        string[] precursorMz = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                        singleSpectrum.precursorMz = double.Parse(precursorMz[1]);
                                    }
                                    catch (Exception e)
                                    {
                                        //singleSpectrum.precursorMz = 0;
                                        throw new MetaMorpheusException("Could not find the mz of precursor : " + e.Message);
                                    }
                                }

                                if (comment[l].Contains("iRT", StringComparison.OrdinalIgnoreCase) || comment[l].Contains("retention time", StringComparison.OrdinalIgnoreCase))
                                {
                                    
                                    try
                                    {
                                        string[] rententionTime = comment[l].Split(new char[] { ':', '=' }).Select(b => b.Trim()).ToArray();
                                        singleSpectrum.rententionTime = double.Parse(rententionTime[1]);
                                    }
                                    catch 
                                    {
                                        //if(e.Message.Contains)
                                        //throw new MetaMorpheusException("Could not find the rentention time : " + e.Message);
                                    }
                                }

                            }
                        }

                        if (lines[j].Contains("peaks", StringComparison.OrdinalIgnoreCase))
                        {
                        var numberOfPeaks = 0;
                        string[] numPeaks = lines[j].Split(":").Select(b => b.Trim()).ToArray();
                        try
                        {
                            numberOfPeaks = int.Parse(numPeaks[1]);
                        }
                        catch
                        {
                            //Console.WriteLine("erroe" + lines[j] + "  " + numPeaks.Length);
                            //foreach (var x in numPeaks)
                            //{
                            //    Console.WriteLine(x);
                            //}
                 
                        }
                            
                            var peaksList = new List<Product>();
                            var peaksWithIntensityList = new Dictionary<Product, double>();
                            for (int k = j + 1; k < j + 1 + numberOfPeaks; k++)
                            {
                                string[] eachPeak = lines[k].Split("\t").Select(b => b.Trim()).ToArray();
                                var experMz = double.Parse(eachPeak[0]);
                                var experIntensity = double.Parse(eachPeak[1]);
                                string[] by = eachPeak[2].Split(new char[] { '/', '\"', 'p' }, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                                var spectrumPeakProductType = by[0].ToCharArray()[0].ToString();
                                var fragmentNumber = (int)Char.GetNumericValue(by[0].ToCharArray()[1]);
            
                                ProductType peakProductType = (ProductType)Enum.Parse(typeof(ProductType), spectrumPeakProductType, true);
                                FragmentationTerminus terminus = (FragmentationTerminus)Enum.Parse(typeof(FragmentationTerminus), "None", true);
                                ProductType.a.ToString();
                                var product = new Product(peakProductType, terminus, experMz, fragmentNumber, 0, 0);
                                //var b = new MatchedFragmentIon(ref product, double.Parse(eachPeak[0]), double.Parse(eachPeak[1]), 1);
                                //b.MassErrorPpm = double.Parse(by[1].Trim());

                                //string[] by = eachPeak[2].Split(new char[] { '/', '\"', 'p'}, StringSplitOptions.RemoveEmptyEntries).Select(b => b.Trim()).ToArray();
                                //b.SpectrumPeakProductType = by[0].ToCharArray()[0].ToString();
                                //b.fragmentNumber = (int)Char.GetNumericValue(by[0].ToCharArray()[1]);
                                ////Console.WriteLine(by[0].ToCharArray()[1]);
                                ////Console.WriteLine(b.fragmentNumber);
                                //b.massErrorPpm = double.Parse(by[1].Trim());
                                peaksList.Add(product);
                                peaksWithIntensityList.Add(product, experIntensity);
                            }
                            singleSpectrum.Peaks = peaksList;
                            singleSpectrum.PeaksWithIntensity = peaksWithIntensityList;
                        }
                    }
                    if (singleSpectrum.Name != null &&  singleSpectrum.Peaks.Count != 0)
                    {
                        SpectralLibraryDictionary.Add(singleSpectrum.Name, singleSpectrum);
                    }
                }
                if (SpectralLibraryDictionary.Count == 0)
                {
                    Console.WriteLine("the library doesn't contain any spectrums!");
                }
         
        }
    }
}
