
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

namespace Test
{
    [TestFixture]
    public static class SpectralLibraryTest
    {
        [Test]
        public static void SpectrumTest()

        {


            string myMzMLTestFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            DbForTask db = new DbForTask(myDatabase, true);
            DbForTask dbOfSpectralLibrary = new DbForTask(spectralLibrary, true, true);
            string spectralLibraryTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralLibrary");
            Directory.CreateDirectory(spectralLibraryTestFolder);
            //run task
            var x = new SearchTask();
            x.RunTask(spectralLibraryTestFolder, new List<DbForTask> { db, dbOfSpectralLibrary }, new List<string> { myMzMLTestFile }, "normal");


            //string myMzMLTestFile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\04-30-13_CAST_Frac5_4uL.raw");
            //string myDatab = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\uniprot-mouse-reviewed-2-5-2018.fasta");
            //string spectralLib = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Mouse_spectralLibrary.msp");
            //DbForTask db = new DbForTask(myDatab, true);
            //DbForTask dbOfSpectralLibrary = new DbForTask(spectralLib, true, true);
            //string spectralLibraryTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralLibrary");
            //Directory.CreateDirectory(spectralLibraryTestFolder);
            //var x = new SearchTask();
            //x.RunTask(spectralLibraryTestFolder, new List<DbForTask> { db, dbOfSpectralLibrary }, new List<string> { myMzMLTestFile1 }, "normal");

            //string pathToExpectedSpectralLibaryResult = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            //string TestSpectralLibraryResultFullPath = Path.Combine(spectralLibraryTestFolder, @"spectralLibrary.msp");

            //var test = y.GetAllScansList();
            //foreach (var ms2scan in test.Where(x => x.MsnOrder != 1))
            //{
            //    var xx = ms2scan.ScanFilter;
            //    Console.WriteLine(xx);
            //}
            //    var _first = test.First();
            //var _second = _first.ScanFilter;
            //Console.WriteLine(_second);
            //Console.WriteLine(_first);

            //Console.WriteLine("y");

            //var expectedFileContents = File.ReadAllLines(pathToExpectedSpectralLibaryResult);
            //var testFileContenets = File.ReadAllLines(TestSpectralLibraryResultFullPath);
            //for (int i = 0; i < expectedFileContents.Length; i++)
            //{
            //    Assert.AreEqual(expectedFileContents[i], testFileContenets[i]);
            //}

            //var y = new LoadAllStaticData(myMzMLTestFile, null, -1);



            //Directory.Delete(spectralLibraryTestFolder, true);
        }
    }
}
