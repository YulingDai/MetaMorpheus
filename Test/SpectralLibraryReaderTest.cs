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
            string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            string spectral= Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            string xyz = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\testofspectralLibrary.msp");
            spectralLibraryReader x = new spectralLibraryReader(xyz);

            //spectralLibraryReader x = new spectralLibraryReader(spectralLibrary);
            foreach(var y in x.SpectralLibraryDictionary)
            {
                //Console.WriteLine(y.Key);
                //Console.WriteLine(y.Value.ToString());
            }

            string myMzMLTestFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            //string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            DbForTask db = new DbForTask(myDatabase, true);
            DbForTask dbOfSpectralLibrary = new DbForTask(spectralLibrary, true, true);
         

            string spectralLibraryTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralLibrary");
            Directory.CreateDirectory(spectralLibraryTestFolder);

            string pathToExpectedSpectralLibaryResult = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            string TestSpectralLibraryResultFullPath = Path.Combine(spectralLibraryTestFolder, @"spectralLibrary.msp");


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
