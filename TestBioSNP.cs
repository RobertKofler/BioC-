/*
"The contents of this file are subject to the Mozilla Public License
Version 1.1 (the "License"); you may not use this file except in
compliance with the License. You may obtain a copy of the License at
http://www.mozilla.org/MPL/

Software distributed under the License is distributed on an "AS IS"
basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
License for the specific language governing rights and limitations
under the License.

The Original Code is the content of this file.

The Initial Developer of the Original Code is Dr. Robert Kofler.
Carrer del Baluard; 73 2°1a
08003 Barcelona
All Rights Reserved.

Contributor(s): Dr. Robert Kofler.

Alternatively, the contents of this file may be used under the terms
of the _____ license (the  "[___] License"), in which case the
provisions of [______] License are applicable instead of those
above. If you wish to allow use of your version of this file only
under the terms of the [____] License and not to allow others to use
your version of this file under the MPL, indicate your decision by
deleting the provisions above and replace them with the notice and
other provisions required by the [___] License. If you do not delete
the provisions above, a recipient may use your version of this file
under either the MPL or the [___] License."
*/
using System;
using System.Collections.Generic;
using System.Text;
using Bio.SNP;
using Bio.Alignment;
using Bio.Seq;

namespace TestBioSNP
{
    using NUnit.Framework;
    class Program
    {
        static void Main(string[] args)
        {
            TestSNP.Start();
        }
    }

    [TestFixture]
    public class TestSNP
    {
        public static void Start()
        {
            TestSNP p = new TestSNP();
            p.Test_SNP();
            p.Test_HalfSNPExtractor();
            p.Test_HalfSNPExtractor_Indel();
            p.Test_FullSNPExtractorSite();
            p.Test_FullSNPExtractorSiteIndel();
            p.Test_QualityAssessementNormal();
            p.Test_QualityAssessement454();
        }

        [Test]
        public void Test_SNP()
        {
            SNP s = new SNP("franz", "fritz", 'T', 'A', 234, 1, 124, 0, false, true, true, 10);
            s.AlleleID = "allele1";
            s.AlignmentQualityTokens = 12;


            Assert.AreEqual(s.DatabaseParent, "franz");
            Assert.AreEqual(s.QueryParent, "fritz");
            Assert.AreEqual(s.DatabaseChar, 'T');
            Assert.AreEqual(s.QueryChar, 'A');
            Assert.AreEqual(s.SiteBasis_Database, 234);
            Assert.AreEqual(s.SiteIndelShift_Database, 1);
            Assert.AreEqual(s.SiteBasis_Query, 124);
            Assert.AreEqual(s.SiteIndelShift_Query, 0);
            Assert.AreEqual(s.SiteID_Database, "234-1");
            Assert.AreEqual(s.SiteID_Query, "124-0");
            Assert.AreEqual(s.Start, 234);
            Assert.AreEqual(s.End, 234);
            Assert.AreEqual(s.Length, 1);
            Assert.AreEqual(s.AlleleID, "allele1");
            Assert.AreEqual(s.Transition, false);
            Assert.AreEqual(s.Transversion, true);
            Assert.AreEqual(s.AlignmentQualityTokens, 12);
            Assert.AreEqual(s.SequenceQualityAtSite, null);
            Assert.AreEqual(s.SequenceQualityInNeighborhood, null);
            Assert.AreEqual(s.PlusPlusStrand, false);
            Assert.AreEqual(s.IsFullSNP, true);
            Assert.AreEqual(s.Valid, true);

            Assert.AreEqual(s.DistanceFromAlignmentEnd, 10);

            s = new SNP("franz", "fritz", 'A', 'A', "1-0", "102-1");
            s.PlusPlusStrand = false;
            s.AlignmentQualityTokens = 3;
            s.DistanceFromAlignmentEnd = 10;
            s.SequenceQualityAtSite = 12;
            s.SequenceQualityInNeighborhood = 20;


            Assert.AreEqual(s.DatabaseParent, "franz");
            Assert.AreEqual(s.QueryParent, "fritz");
            Assert.AreEqual(s.DatabaseChar, 'A');
            Assert.AreEqual(s.QueryChar, 'A');
            Assert.AreEqual(s.SiteBasis_Database, 1);
            Assert.AreEqual(s.SiteIndelShift_Database, 0);
            Assert.AreEqual(s.SiteID_Database, "1-0");
            Assert.AreEqual(s.SiteBasis_Query, 102);
            Assert.AreEqual(s.SiteIndelShift_Query, 1);
            Assert.AreEqual(s.SiteID_Query, "102-1");
            Assert.AreEqual(s.Start, 1);
            Assert.AreEqual(s.End, 1);
            Assert.AreEqual(s.Length, 1);
            Assert.AreEqual(s.AlleleID, null);
            Assert.AreEqual(s.Transition, false);
            Assert.AreEqual(s.Transversion, false);
            Assert.AreEqual(s.PlusPlusStrand, false);
            Assert.AreEqual(s.AlignmentQualityTokens, 3);
            Assert.AreEqual(s.SequenceQualityAtSite, 12);
            Assert.AreEqual(s.SequenceQualityInNeighborhood, 20);
            Assert.AreEqual(s.DistanceFromAlignmentEnd, 10);

        }

        [Test]
        public void Test_HalfSNPExtractor()
        {

            ISNPExtractor extractor = new HalfSNPExtractor();
            PairwiseNucleotideSequenceAlignment pw;
            List<SNP> results;

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "102-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");
            Assert.AreEqual(results[0].IsFullSNP, false);
            Assert.AreEqual(results[0].Valid, true);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAATC"), new Sequence_ByteArray("CAAACG")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'C');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'T');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "104-0");
            Assert.AreEqual(results[1].SiteID_Query, "14-0");
            Assert.AreEqual(results[2].DatabaseChar, 'C');
            Assert.AreEqual(results[2].QueryChar, 'G');
            Assert.AreEqual(results[2].IsFullSNP, false);
            Assert.AreEqual(results[2].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 1);
            Assert.AreEqual(results[2].DistanceFromAlignmentEnd, 0);


            extractor = new HalfSNPExtractor();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'C');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CAAACAA"), new Sequence_ByteArray("GAAAGAA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);

            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'G');
            Assert.AreEqual(results[1].SiteID_Database, "104-0");
            Assert.AreEqual(results[1].SiteID_Query, "14-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 2);


            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);

            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 2);


            ///Test gaps
            ///
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAGGG"), new Sequence_ByteArray("AAA-CA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            extractor = new HalfSNPExtractor();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'G');
            Assert.AreEqual(results[1].QueryChar, 'A');
            Assert.AreEqual(results[1].SiteID_Database, "105-0");
            Assert.AreEqual(results[1].SiteID_Query, "14-0");

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CAAAAA"), new Sequence_ByteArray("GC-AAA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAA--CAT"), new Sequence_ByteArray("A-AAACAC")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'T');
            Assert.AreEqual(results[0].QueryChar, 'C');
            Assert.AreEqual(results[0].SiteID_Database, "105-0");
            Assert.AreEqual(results[0].SiteID_Query, "16-0");

            ///Round 2 Composite Alignments
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAAAC"), new Sequence_ByteArray("AAA-CG")), "fritz", "franz", 200, 20, 205, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            results = extractor.GetSNPs(com);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "102-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'A');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "204-0");
            Assert.AreEqual(results[1].SiteID_Query, "23-0");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 1);

            Assert.AreEqual(results[2].DatabaseParent, "fritz");
            Assert.AreEqual(results[2].QueryParent, "franz");
            Assert.AreEqual(results[2].DatabaseChar, 'C');
            Assert.AreEqual(results[2].QueryChar, 'G');
            Assert.AreEqual(results[2].SiteID_Database, "205-0");
            Assert.AreEqual(results[2].SiteID_Query, "24-0");
            Assert.AreEqual(results[2].DistanceFromAlignmentEnd, 0);

        }

        [Test]
        public void Test_HalfSNPExtractor_Indel()
        {

            ISNPExtractor extractor = new HalfSNPExtractor_Indel();
            PairwiseNucleotideSequenceAlignment pw;
            List<SNP> results;

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "102-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");
            Assert.AreEqual(results[0].IsFullSNP, false);
            Assert.AreEqual(results[0].Valid, true);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAATC"), new Sequence_ByteArray("CAAACG")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'C');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'T');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "104-0");
            Assert.AreEqual(results[1].SiteID_Query, "14-0");
            Assert.AreEqual(results[2].DatabaseChar, 'C');
            Assert.AreEqual(results[2].QueryChar, 'G');
            Assert.AreEqual(results[2].IsFullSNP, false);
            Assert.AreEqual(results[2].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 1);
            Assert.AreEqual(results[2].DistanceFromAlignmentEnd, 0);


            extractor = new HalfSNPExtractor_Indel();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'C');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CAAACAA"), new Sequence_ByteArray("GAAAGAA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);

            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'G');
            Assert.AreEqual(results[1].SiteID_Database, "104-0");
            Assert.AreEqual(results[1].SiteID_Query, "14-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 2);


            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);

            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 2);


            ///Test gaps
            ///
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAGGG"), new Sequence_ByteArray("AAA-CA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            extractor = new HalfSNPExtractor_Indel();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'G');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "103-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-1");

            Assert.AreEqual(results[2].DatabaseParent, "fritz");
            Assert.AreEqual(results[2].QueryParent, "franz");
            Assert.AreEqual(results[2].DatabaseChar, 'G');
            Assert.AreEqual(results[2].QueryChar, 'A');
            Assert.AreEqual(results[2].SiteID_Database, "105-0");
            Assert.AreEqual(results[2].SiteID_Query, "14-0");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAGGG"), new Sequence_ByteArray("AAT--G")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            extractor = new HalfSNPExtractor_Indel();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'T');
            Assert.AreEqual(results[0].SiteID_Database, "102-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'G');
            Assert.AreEqual(results[1].QueryChar, '-');
            Assert.AreEqual(results[1].SiteID_Database, "103-0");
            Assert.AreEqual(results[1].SiteID_Query, "12-1");

            Assert.AreEqual(results[2].DatabaseParent, "fritz");
            Assert.AreEqual(results[2].QueryParent, "franz");
            Assert.AreEqual(results[2].DatabaseChar, 'G');
            Assert.AreEqual(results[2].QueryChar, '-');
            Assert.AreEqual(results[2].SiteID_Database, "104-0");
            Assert.AreEqual(results[2].SiteID_Query, "12-2");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AA---G"), new Sequence_ByteArray("AAAGG-")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            extractor = new HalfSNPExtractor_Indel();
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 4);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, '-');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-1");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, '-');
            Assert.AreEqual(results[1].QueryChar, 'G');
            Assert.AreEqual(results[1].SiteID_Database, "101-2");
            Assert.AreEqual(results[1].SiteID_Query, "13-0");

            Assert.AreEqual(results[2].DatabaseParent, "fritz");
            Assert.AreEqual(results[2].QueryParent, "franz");
            Assert.AreEqual(results[2].DatabaseChar, '-');
            Assert.AreEqual(results[2].QueryChar, 'G');
            Assert.AreEqual(results[2].SiteID_Database, "101-3");
            Assert.AreEqual(results[2].SiteID_Query, "14-0");

            Assert.AreEqual(results[3].DatabaseParent, "fritz");
            Assert.AreEqual(results[3].QueryParent, "franz");
            Assert.AreEqual(results[3].DatabaseChar, 'G');
            Assert.AreEqual(results[3].QueryChar, '-');
            Assert.AreEqual(results[3].SiteID_Database, "102-0");
            Assert.AreEqual(results[3].SiteID_Query, "14-1");

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CAAAAA"), new Sequence_ByteArray("GC-AAA")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "100-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAA--CAT"), new Sequence_ByteArray("A-AAACAC")), "fritz", "franz", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 4);
            Assert.AreEqual(results[3].DatabaseParent, "fritz");
            Assert.AreEqual(results[3].QueryParent, "franz");
            Assert.AreEqual(results[3].DatabaseChar, 'T');
            Assert.AreEqual(results[3].QueryChar, 'C');
            Assert.AreEqual(results[3].SiteID_Database, "105-0");
            Assert.AreEqual(results[3].SiteID_Query, "16-0");

            ///Round 2 Composite Alignments
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAAAC"), new Sequence_ByteArray("AAA-CG")), "fritz", "franz", 200, 20, 205, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            results = extractor.GetSNPs(com);
            Assert.AreEqual(results.Count, 4);
            Assert.AreEqual(results[0].DatabaseParent, "fritz");
            Assert.AreEqual(results[0].QueryParent, "franz");
            Assert.AreEqual(results[0].DatabaseChar, 'C');
            Assert.AreEqual(results[0].QueryChar, 'G');
            Assert.AreEqual(results[0].SiteID_Database, "102-0");
            Assert.AreEqual(results[0].SiteID_Query, "12-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);

            Assert.AreEqual(results[1].DatabaseParent, "fritz");
            Assert.AreEqual(results[1].QueryParent, "franz");
            Assert.AreEqual(results[1].DatabaseChar, 'A');
            Assert.AreEqual(results[1].QueryChar, '-');
            Assert.AreEqual(results[1].SiteID_Database, "203-0");
            Assert.AreEqual(results[1].SiteID_Query, "22-1");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 2);

            Assert.AreEqual(results[2].DatabaseParent, "fritz");
            Assert.AreEqual(results[2].QueryParent, "franz");
            Assert.AreEqual(results[2].DatabaseChar, 'A');
            Assert.AreEqual(results[2].QueryChar, 'C');
            Assert.AreEqual(results[2].SiteID_Database, "204-0");
            Assert.AreEqual(results[2].SiteID_Query, "23-0");
            Assert.AreEqual(results[2].DistanceFromAlignmentEnd, 1);

            Assert.AreEqual(results[3].DatabaseParent, "fritz");
            Assert.AreEqual(results[3].QueryParent, "franz");
            Assert.AreEqual(results[3].DatabaseChar, 'C');
            Assert.AreEqual(results[3].QueryChar, 'G');
            Assert.AreEqual(results[3].SiteID_Database, "205-0");
            Assert.AreEqual(results[3].SiteID_Query, "24-0");
            Assert.AreEqual(results[3].DistanceFromAlignmentEnd, 0);



        }

        [Test]
        public void Test_FullSNPExtractorSite()
        {

            List<SNPSite> sites = new List<SNPSite>();
            sites.Add(new SNPSite("franz", "99-0"));
            sites.Add(new SNPSite("franz", "101-0"));
            sites.Add(new SNPSite("franz", "112-0"));

            Dictionary<string, List<SNPSite>> dict = new Dictionary<string, List<SNPSite>>();
            dict.Add("franz", sites);

            ISNPExtractor extractor = new FullSNPExtractorSite(dict);
            PairwiseNucleotideSequenceAlignment pw;
            List<SNP> results;

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "11-0");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AACTTTTTAAAAC"), new Sequence_ByteArray("AAGTTTTTAAAAC")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "11-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdinand");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "112-0");
            Assert.AreEqual(results[1].SiteID_Query, "22-0");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC----TTTTTAAAAC"), new Sequence_ByteArray("A--AAAGTTTTTAAAAC")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdinand");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "112-0");
            Assert.AreEqual(results[1].SiteID_Query, "24-0");
            Assert.AreEqual(results[1].IsFullSNP, true);
            Assert.AreEqual(results[1].Valid, true);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);

            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-1");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, false);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AYG")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, false);
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'Y');
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("A-G")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-1");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, false);

            //Round2

            sites = new List<SNPSite>();
            sites.Add(new SNPSite("franz", "99-0"));
            sites.Add(new SNPSite("franz", "101-2"));
            sites.Add(new SNPSite("franz", "112-0"));
            dict = new Dictionary<string, List<SNPSite>>();
            dict.Add("franz", sites);
            extractor = new FullSNPExtractorSite(dict);


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdinand", 99, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "99-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);

            ///Round 3 Composite Alignments
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdl", 99, 10, 101, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAAAC"), new Sequence_ByteArray("AAA-CG")), "franz", "ferdl", 107, 20, 112, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            results = extractor.GetSNPs(com);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdl");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "99-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);



            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdl");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'G');
            Assert.AreEqual(results[1].SiteID_Database, "112-0");
            Assert.AreEqual(results[1].SiteID_Query, "24-0");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);
        }

        [Test]
        public void Test_FullSNPExtractorSiteIndel()
        {

            List<SNPSite> sites = new List<SNPSite>();
            sites.Add(new SNPSite("franz", "99-0"));
            sites.Add(new SNPSite("franz", "101-0"));
            sites.Add(new SNPSite("franz", "112-0"));

            Dictionary<string, List<SNPSite>> dict = new Dictionary<string, List<SNPSite>>();
            dict.Add("franz", sites);

            ISNPExtractor extractor = new FullSNPExtractorSite_Indel(dict);
            PairwiseNucleotideSequenceAlignment pw;
            List<SNP> results;

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "11-0");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AACTTTTTAAAAC"), new Sequence_ByteArray("AAGTTTTTAAAAC")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "11-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdinand");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "112-0");
            Assert.AreEqual(results[1].SiteID_Query, "22-0");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC----TTTTTAAAAC"), new Sequence_ByteArray("A--AAAGTTTTTAAAAC")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = false;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdinand");
            Assert.AreEqual(results[1].DatabaseChar, 'C');
            Assert.AreEqual(results[1].QueryChar, 'C');
            Assert.AreEqual(results[1].SiteID_Database, "112-0");
            Assert.AreEqual(results[1].SiteID_Query, "24-0");
            Assert.AreEqual(results[1].IsFullSNP, true);
            Assert.AreEqual(results[1].Valid, true);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);

            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-1");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AYG")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, false);
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'Y');
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 1);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("A-G")), "franz", "ferdinand", 100, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-1");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);

            //Round2

            sites = new List<SNPSite>();
            sites.Add(new SNPSite("franz", "99-0"));
            sites.Add(new SNPSite("franz", "101-2"));
            sites.Add(new SNPSite("franz", "112-0"));
            dict = new Dictionary<string, List<SNPSite>>();
            dict.Add("franz", sites);
            extractor = new FullSNPExtractorSite_Indel(dict);


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdinand", 99, 10, 150, 50);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "99-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);

            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdinand");
            Assert.AreEqual(results[1].DatabaseChar, '-');
            Assert.AreEqual(results[1].QueryChar, '-');
            Assert.AreEqual(results[1].SiteID_Database, "101-2");
            Assert.AreEqual(results[1].SiteID_Query, "12-2");
            Assert.AreEqual(results[1].IsFullSNP, true);
            Assert.AreEqual(results[1].Valid, true);
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);

            ///Round 3 Composite Alignments
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "franz", "ferdl", 99, 10, 101, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAAAC"), new Sequence_ByteArray("AAA-CG")), "franz", "ferdl", 107, 20, 112, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            results = extractor.GetSNPs(com);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdl");
            Assert.AreEqual(results[0].DatabaseChar, 'A');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "99-0");
            Assert.AreEqual(results[0].SiteID_Query, "10-0");
            Assert.AreEqual(results[0].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);

            Assert.AreEqual(results[1].DatabaseParent, "franz");
            Assert.AreEqual(results[1].QueryParent, "ferdl");
            Assert.AreEqual(results[1].DatabaseChar, '-');
            Assert.AreEqual(results[1].QueryChar, '-');
            Assert.AreEqual(results[1].SiteID_Database, "101-2");
            Assert.AreEqual(results[1].SiteID_Query, "12-2");
            Assert.AreEqual(results[1].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[1].IsFullSNP, true);
            Assert.AreEqual(results[1].Valid, true);

            Assert.AreEqual(results[2].DatabaseParent, "franz");
            Assert.AreEqual(results[2].QueryParent, "ferdl");
            Assert.AreEqual(results[2].DatabaseChar, 'C');
            Assert.AreEqual(results[2].QueryChar, 'G');
            Assert.AreEqual(results[2].SiteID_Database, "112-0");
            Assert.AreEqual(results[2].SiteID_Query, "24-0");
            Assert.AreEqual(results[2].DistanceFromAlignmentEnd, 0);
            Assert.AreEqual(results[2].IsFullSNP, true);
            Assert.AreEqual(results[2].Valid, true);



            //Round 4
            sites = new List<SNPSite>();
            sites.Add(new SNPSite("franz", "101-4"));
            dict = new Dictionary<string, List<SNPSite>>();
            dict.Add("franz", sites);
            extractor = new FullSNPExtractorSite_Indel(dict);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AA----A"), new Sequence_ByteArray("AAAAAAA")), "franz", "ferdinand", 100, 10, 106, 16);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, '-');
            Assert.AreEqual(results[0].QueryChar, 'A');
            Assert.AreEqual(results[0].SiteID_Database, "101-4");
            Assert.AreEqual(results[0].SiteID_Query, "15-0");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AA--A"), new Sequence_ByteArray("AAAAA")), "franz", "ferdinand", 100, 10, 106, 16);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, '-');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-4");
            Assert.AreEqual(results[0].SiteID_Query, "13-2");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AA---A"), new Sequence_ByteArray("AA-A-A")), "franz", "ferdinand", 100, 10, 106, 16);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, '-');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-4");
            Assert.AreEqual(results[0].SiteID_Query, "12-1");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAA"), new Sequence_ByteArray("AAA")), "franz", "ferdinand", 101, 10, 103, 13);
            pw.PlusPlusStrand = true;
            results = extractor.GetSNPs(pw);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseParent, "franz");
            Assert.AreEqual(results[0].QueryParent, "ferdinand");
            Assert.AreEqual(results[0].DatabaseChar, '-');
            Assert.AreEqual(results[0].QueryChar, '-');
            Assert.AreEqual(results[0].SiteID_Database, "101-4");
            Assert.AreEqual(results[0].SiteID_Query, "10-4");
            Assert.AreEqual(results[0].IsFullSNP, true);
            Assert.AreEqual(results[0].Valid, true);
        }

        [Test]
        public void Test_QualityAssessementNormal()
        {

            //Round1
            ISNPQualityAlignment quality;
            PairwiseNucleotideSequenceAlignment pw;
            SNP s;

            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 21, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            int qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 1);

            s = new SNP("fritz", "franz", 'A', 'C', 22, 0, 12, 0, true, true, true, 0);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 0);

            quality = new SNPAlignmentQualityAssessment(3);
            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 1);

            quality = new SNPAlignmentQualityAssessment(3);
            s = new SNP("fritz", "franz", 'A', 'C', 21, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 0);

            quality = new SNPAlignmentQualityAssessment(3);
            s = new SNP("fritz", "franz", 'A', 'C', 18, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 0);

            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 18, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 1);

            quality = new SNPAlignmentQualityAssessment(7);
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA---ATTTTTT"), new Sequence_ByteArray("AAAACCCATTTTTT")), "fritz", "franz", 200, 20, 210, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            s = new SNP("fritz", "franz", 'A', 'C', 204, 0, 12, 0, true, true, true, 0);
            qual = quality.GetLowQualityToken(s, com);
            Assert.AreEqual(qual, 3);

            temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA--ATTTTTT"), new Sequence_ByteArray("AAAACCATTTTTT")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA"), new Sequence_ByteArray("AAAA")), "fritz", "franz", 200, 20, 210, 25, true));
            com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            s = new SNP("fritz", "franz", 'A', 'C', 104, 0, 12, 0, true, true, true, 0);
            qual = quality.GetLowQualityToken(s, com);
            Assert.AreEqual(qual, 2);



            //Round 2
            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 18, 1, 12, 1, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TTAA"), new Sequence_ByteArray("AAATTAA")), "fritz", "franz", 18, 10, 23, 15);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 1);

            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 19, 3, 12, 1, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC---AAAA"), new Sequence_ByteArray("ACAAAAAAA")), "fritz", "franz", 18, 10, 23, 18);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 3);

            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 19, 2, 12, 1, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC---AAAA"), new Sequence_ByteArray("A-AAAAAAA")), "fritz", "franz", 18, 10, 23, 18);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 4);


            quality = new SNPAlignmentQualityAssessment(5);
            s = new SNP("fritz", "franz", 'A', 'C', 19, 4, 12, 1, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC---AAAA"), new Sequence_ByteArray("ACAAAAAAA")), "fritz", "franz", 18, 10, 23, 18);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 3);


            //Round 3

        }

        [Test]
        public void Test_QualityAssessement454()
        {
            ISNPQualityAlignment quality;
            PairwiseNucleotideSequenceAlignment pw;
            SNP s;

            quality = new SNPAlignmentQualityAssessment_454(5);


            s = new SNP("fritz", "franz", 'A', 'C', 26, 0, 18, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AGAGCCAATAACCAGAG"), new Sequence_ByteArray("AGAGCCAATAACCAGAG")), "fritz", "franz", 18, 10, 34, 26, true);
            int qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 8);

            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ACTTAA"), new Sequence_ByteArray("AATTAA")), "fritz", "franz", 18, 10, 23, 15, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 3);

            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ACCTTAA"), new Sequence_ByteArray("AA-TTAA")), "fritz", "franz", 18, 10, 23, 15, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 4);

            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAACCTTAA"), new Sequence_ByteArray("AAAAA-TTAA")), "fritz", "franz", 15, 10, 23, 15, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 7);


            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TAAAA"), new Sequence_ByteArray("AACTAAAA")), "fritz", "franz", 18, 10, 23, 15, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 5);

            s = new SNP("fritz", "franz", 'A', 'C', 20, 0, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC-TAAAA"), new Sequence_ByteArray("A-CTAAAA")), "fritz", "franz", 18, 10, 23, 15, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 5);


            quality = new SNPAlignmentQualityAssessment_454(7);
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAC"), new Sequence_ByteArray("AAG")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA---ATTTTTT"), new Sequence_ByteArray("AAAACCCATTTTTT")), "fritz", "franz", 200, 20, 210, 25, true));
            CompositePairwiseNucleotideSequenceAlignment com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            s = new SNP("fritz", "franz", 'A', 'C', 204, 0, 12, 0, true, true, true, 0);
            qual = quality.GetLowQualityToken(s, com);
            Assert.AreEqual(qual, 14);

            temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA--ATTTTTT"), new Sequence_ByteArray("AAAACCATTTTTT")), "fritz", "franz", 100, 10, 110, 13, true));
            temp.Add(new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAA"), new Sequence_ByteArray("AAAA")), "fritz", "franz", 200, 20, 210, 25, true));
            com = new CompositePairwiseNucleotideSequenceAlignment(temp);
            s = new SNP("fritz", "franz", 'A', 'C', 104, 0, 12, 0, true, true, true, 0);
            qual = quality.GetLowQualityToken(s, com);
            Assert.AreEqual(qual, 13);



            //Round 2
            quality = new SNPAlignmentQualityAssessment_454(5);
            s = new SNP("fritz", "franz", 'A', 'C', 19, 1, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC----AATTTT"), new Sequence_ByteArray("ACTCTCAATTTT")), "fritz", "franz", 18, 10, 25, 20, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 3);

            s = new SNP("fritz", "franz", 'A', 'C', 19, 2, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC----AATTTT"), new Sequence_ByteArray("ACTCTCAATTTT")), "fritz", "franz", 18, 10, 25, 20, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 5);

            s = new SNP("fritz", "franz", 'A', 'C', 19, 3, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC----AATTTT"), new Sequence_ByteArray("ACTCTCAATTTT")), "fritz", "franz", 18, 10, 25, 20, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 5);

            s = new SNP("fritz", "franz", 'A', 'C', 19, 4, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC----AATTTT"), new Sequence_ByteArray("ACTCTCAATTTT")), "fritz", "franz", 18, 10, 25, 20, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 7);

            s = new SNP("fritz", "franz", 'A', 'C', 19, 5, 12, 0, true, true, true, 0);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AC----AATTTT"), new Sequence_ByteArray("ACTCTCAATTTT")), "fritz", "franz", 18, 10, 25, 20, true);
            qual = quality.GetLowQualityToken(s, pw);
            Assert.AreEqual(qual, 7);




        }

    }
}
