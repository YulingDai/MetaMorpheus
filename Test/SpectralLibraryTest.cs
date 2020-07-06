
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
        //// helper method to compare  files
        //public static string GetFileHash(string filename)
        //{
        //    var hash = new SHA256Managed();
        //    var clearBytes = File.ReadAllBytes(filename);
        //    var hashedBytes = hash.ComputeHash(clearBytes);
        //    return ConvertBytesToHex(hashedBytes);
        //}

        //// helper method to compare 2 files
        //public static string ConvertBytesToHex(byte[] bytes)
        //{
        //    var sb = new StringBuilder();

        //    for (var i = 0; i < bytes.Length; i++)
        //    {
        //        sb.Append(bytes[i].ToString("x"));
        //    }
        //    return sb.ToString();
        //}


        [Test]
        public static void SpectrumTest()
        {
            string myMzMLTestFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string spectralLibrary = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            DbForTask db = new DbForTask(myDatabase, true);

            string spectralLibraryTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralLibrary");
            Directory.CreateDirectory(spectralLibraryTestFolder);

            string pathToExpectedSpectralLibaryResult = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\spectralLibrary.msp");
            string TestSpectralLibraryResultFullPath = Path.Combine(spectralLibraryTestFolder, @"spectralLibrary.msp");

            
            //run task
            var x = new SearchTask();
            
   
            x.RunTask(spectralLibraryTestFolder, new List<DbForTask> { db }, new List<string> { myMzMLTestFile }, spectralLibrary, "normal");
            //x.RunTask(spectralLibraryTestFolder, new List<DbForTask> { db }, new List<string> { myMzMLTestFile }, spectralLibrary, "normal");
            //new SearchTask().RunTask(spectralLibraryTestFolder, new List<DbForTask> { db }, new List<string> { myMzMLTestFile }, "normal");
            // ThermoRawFileReaderData y = ThermoRawFileReaderData.LoadAllStaticData(myMzMLTestFile,null,-1);
            //Console.WriteLine(y.GetAllScansList());


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
