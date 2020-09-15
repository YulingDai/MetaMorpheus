using System.Collections.Generic;
using System.Text;
using TaskLayer;
using NUnit.Framework;
using System.IO;
using System.Security.Cryptography;
using System;
using ThermoRawFileReader;
using System.Linq;
using MassSpectrometry;
using ThermoFisher.CommonCore.Data.Interfaces;
using EngineLayer.spectralLibrarySearch;

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryReaderTest
    {
        [Test]
        public static void SpectralReaderTest()
        {
            string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\HelaMannSpectralLibrary.msp");
            string spectral= Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\newspectralLibrary.msp");
            string xyz = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\testofspectralLibrary.msp");
            //spectralLibraryReader x = new spectralLibraryReader(spectral);
            //  var x = new int[] { 15, 20, 11, 30, 33, 22, 37, 20, 29, 35, 8, 10, 15, 25 };
            //var y = new double[] { 10, 21, 5, 7, 3, 8, 10, 14, 10 };
            //Console.WriteLine(x.Sum() / x.Length);
            //Console.WriteLine(y.Sum()/y.Length);
            //int a = 0;
            //foreach(var z in y)
            //{
            //    a = a + z * z;
            //}
            //Console.WriteLine(a);
            //Console.WriteLine(y.Sum()* y.Sum());
            //var x = new int[] {1870,1324,1446,1325,1759,1652,1364, 1515, 1065};
            //var y = new int[] {1121,408, 184, 16, 741, 170, 991, 711, 734, 202, 893, 742, 335, 444};

            //var meanx = x.Sum() / x.Length;
            //var meany = y.Sum() / x.Length;
            //Console.WriteLine(x.Sum() / x.Length);
            //Console.WriteLine(y.Sum() / y.Length);
            //int b = 0;
            //foreach (var z in x)
            //{
            //    b = b + (z-meanx) * (z - meanx);
            //}
            //double a = 0;
            //foreach (var z in y)
            //{
            //    a = a + (z-meany) * (z - meany);
            //}
            //Console.WriteLine((a/(y.Length-1)));
            //Console.WriteLine(y.Sum() * y.Sum());

            String a = "EVEDPQVEQLELGGSPGDLQTLALEVARQ";
            String b = "LPVNSPMTKGDTKVMKCVLEVISDSLSKPSPMPVSPECLETLQGDERILSILRHQNLLKELQDLALQGAKERAQQPLKQQQPPKQQQQQQQQQQQEQQHSSFE" +
                "DELSEVFENQSPDAKHRDAAAEVPSRDTMEKRKDSDKGQQDGFEATTEGPRPQAFPEPNQESPMMGDSESPGEDTATNTQSPTSLPSQEHVDPQATGDSERGLSAQQQARK" +
                "AKQEEKEEEEEEEAVAREKAGPEEVPTAASSSHFHAGYKAIQKDDGQSDSQAVDGDGKTEASEALPSEGKGELEHSQQEEDGEEAMVGTPQGLFPQGGKGRELEHKQEEEEE" +
                "EEERLSREWEDKRWSRMDQLAKELTAEKRLEGEDDPDRSMKLSFRTRAYGFRDPGPQLRRGWRPSSREDSVEARSDFEEKKEEEGSANRRAEDQELESLSAIEAELEKVAHQLQALRRG";
            List<String> insulin1_C_peptide = new List<string>();
            List<String> Chromogranin_A = new List<string>();
            List<String> Combination_insulin1_C_peptide_Chromogranin_A = new List<string>();
            foreach (char c in a)
            {
                insulin1_C_peptide.Add(c.ToString());
            }
            foreach (char c in b)
            {
                Chromogranin_A.Add(c.ToString());
            }
            Console.WriteLine(insulin1_C_peptide.Count + "  " + Chromogranin_A.Count);
            for(int i = 0; i < insulin1_C_peptide.Count; i++)
            {
                for(int j = 0; j < Chromogranin_A.Count; j++)
                {
                    Combination_insulin1_C_peptide_Chromogranin_A.Add(insulin1_C_peptide[i] + Chromogranin_A[j]);
                }
            }
            Console.WriteLine(Combination_insulin1_C_peptide_Chromogranin_A.Count);
            foreach (var p in Combination_insulin1_C_peptide_Chromogranin_A)
            {
                Console.WriteLine(p);
            }
            //int i = 0;
            //int j = 0;

            //while(i!= insulin1_C_peptide.Count && j!= Chromogranin_A.Count)
            //{

            //}





            //spectralLibraryReader x = new spectralLibraryReader(spectralLibrary);
            //Console.WriteLine(x.SpectralLibraryDictionary.Count);
            //foreach (var y in x.SpectralLibraryDictionary)
            //{
            //    Console.WriteLine(y.Key);
            //   foreach(var z in y.Value.PeaksWithIntensity)
            //    {
            //        Console.WriteLine(z.Key.NeutralMass + "   " + z.Value + "   " + z.Key.ProductType.ToString() + z.Key.FragmentNumber);
            //    }



            //}

            //string myMzMLTestFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            //string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            ////string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            //DbForTask db = new DbForTask(myDatabase, true);
            //DbForTask dbOfSpectralLibrary = new DbForTask(spectralLibrary, true, true);


            //string spectralLibraryTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralLibrary");
            //Directory.CreateDirectory(spectralLibraryTestFolder);

            //string pathToExpectedSpectralLibaryResult = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            //string TestSpectralLibraryResultFullPath = Path.Combine(spectralLibraryTestFolder, @"spectralLibrary.msp");


            //run task
            //var x = new SearchTask();
            //Console.WriteLine("0000");
            //x.RunTask(spectralLibraryTestFolder, new List<DbForTask> { db, dbOfSpectralLibrary }, new List<string> { myMzMLTestFile }, "normal");
            //spectralLibraryReader y = new spectralLibraryReader(spectral);
            //var search = new ClassicSearchOfSpectralLibrary(x.spectrums, y.spectrums, 0.05, 0.02, 5);
            //Console.WriteLine(search.SpectralLibrarySearchResults.Length);
            //foreach (SpectralLibrarySearchResults z in search.SpectralLibrarySearchResults)
            //{
            //    Console.WriteLine(z.TheExperimentalSpectrum.ToString() + "matches:" );
            //    foreach(SpectralLibrarayMatch A in z.SpectralLibrarayMatchs)
            //    {
            //        Console.WriteLine("search results: ");
            //        Console.WriteLine(A.MatchedSpectrumFromLibrary.Name);
            //        Console.WriteLine(A.MatchedSpectrumFromLibrary.precursorMz);
            //        Console.WriteLine(A.MatchScore);
            //    }

            //}
            //Assert.AreEqual(y.spectrums[0].Name, "KAPAGGAADAAAK");
            //Assert.AreEqual(search.SpectralLibrarySearchResults[0].TheExperimentalSpectrum.Name, "KAPAGGAADAAAK");
            //Assert.AreEqual(search.SpectralLibrarySearchResults[0].SpectralLibrarayMatchs[0].MatchedSpectrumFromLibrary.Name, "KAPAGGAADAAAK");


        }

       
    }
}
