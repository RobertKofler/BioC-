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
using System.Diagnostics;
using System.IO;

namespace TestBioAln
{

    using Bio.Alignment;
    using Bio.Alignment.Misc;
    using Bio.Blast;
    using NUnit.Framework;
    using Bio.Seq;
    using Bio.Seq.IO;
    using Bio.Alignment.IO;
    using Bio.Blast.Misc;
    

    public class Program
    {
        public static void Main()
        {
            Aln_Basic.Start();
            Aln_Misc.Start();
            Aln_IO.Start();
            Aln_Blast.Start();
        }
    }

    [TestFixture]
    public class Aln_Basic
    {
        public static void Start()
        {
            Aln_Basic p = new Aln_Basic();
            p.Test_PairwiseAlignmentFormater_Blast();
            p.Test_PairwiseAlignmentFormater_WaterEmboss();
            p.Test_PairwiseAlignment();
            p.Test_PairwiseAlignmentBuilder();
            p.Test_SmithWatermanGotoh();
            p.Test_NeedlemanWunschGotoh();
            p.Test_SmithWatermanGotoh_454P();
            p.Test_SmithWatermanGotoh_454M();
            p.Test_PairwiseNucleotideSequenceAlignment();

        }

        [Test]
        public void Test_PairwiseAlignmentFormater_Blast()
        {
            IAlignmentFormater af = new PairwiseAlignmentFormater_Blast();
            ISequenceContainer simSeq= af.GetSimilaritySequence(new PairwiseAlignment(new Sequence_ByteArray("AATTAATTCTTT"), new Sequence_ByteArray("AATTAATT-TTT")));
            Assert.AreEqual(simSeq.ToString(), "|||||||| |||");

            simSeq = af.GetSimilaritySequence(new PairwiseAlignment(new Sequence_ByteArray("AXTTAATTCTTT"), new Sequence_ByteArray("AATTAATT-TNT")));
            Assert.AreEqual(simSeq.ToString(), "| |||||| | |");



        }

        [Test]
        public void Test_PairwiseAlignmentFormater_WaterEmboss()
        {
            IAlignmentFormater af = new PairwiseAlignmentFormater_WaterEmboss();
            ISequenceContainer simSeq = af.GetSimilaritySequence(new PairwiseAlignment(new Sequence_ByteArray("AATTAATTCTTT"), new Sequence_ByteArray("AATTAATT-TTT")));
            Assert.AreEqual(simSeq.ToString(), "|||||||| |||");

            simSeq = af.GetSimilaritySequence(new PairwiseAlignment(new Sequence_ByteArray("AXTTAATTCTTT"), new Sequence_ByteArray("AATTAATT-TNT")));
            Assert.AreEqual(simSeq.ToString(), "| |||||| | |");

            simSeq = af.GetSimilaritySequence(new PairwiseAlignment(new Sequence_ByteArray("CXTTAATTCTTT"), new Sequence_ByteArray("AATTAATT-CNT")));
            Assert.AreEqual(simSeq.ToString(), ". |||||| . |");



        }

        [Test]
        public void Test_PairwiseAlignment()
        {
            PairwiseAlignment alignment = new PairwiseAlignment(new Sequence_ByteArray("AATTAATTCTTT"), new Sequence_ByteArray("AATTAATT-TTT"));

            Assert.AreEqual(alignment.Length, 12);
            Assert.AreEqual(alignment.SimilaritySequence.ToString(), "|||||||| |||");
            Assert.AreEqual(alignment.DatabaseSequence.ToString(), "AATTAATTCTTT");
            Assert.AreEqual(alignment.QuerySequence.ToString(), "AATTAATT-TTT");
            Assert.AreEqual(alignment.QuerySequence_WithoutIndels.ToString(), "AATTAATTTTT");
            Assert.AreEqual(alignment[1], new char[] { 'A', '|', 'A' });

            alignment = alignment.GetSubAlignment(9, 3);
            Assert.AreEqual(alignment.Length, 3);
            Assert.AreEqual(alignment.DatabaseSequence.ToString(), "TTT");
            Assert.AreEqual(alignment.QuerySequence.ToString(), "TTT");
            Assert.AreEqual(alignment.SimilaritySequence.ToString(), "|||");

            alignment = new PairwiseAlignment(new Sequence_ByteArray("AATTA---ATTCCCC"), new Sequence_ByteArray("AATTAATT---CCCC"));
            Assert.AreEqual(alignment.DatabaseSequence_WithoutIndels.ToString(), "AATTAATTCCCC");
            Assert.AreEqual(alignment.DatabaseSequence.ToString(), "AATTA---ATTCCCC");
            Assert.AreEqual(alignment.QuerySequence.ToString(), "AATTAATT---CCCC");
            Assert.AreEqual(alignment.SimilaritySequence.ToString(), "|||||      ||||");
            Assert.AreEqual(alignment.QuerySequence_WithoutIndels.ToString(), "AATTAATTCCCC");
            Assert.AreEqual(alignment.ToString(), "AATTA---ATTCCCC\r\n|||||      ||||\r\nAATTAATT---CCCC\r\n");

            alignment = new PairwiseAlignment(new Sequence_ByteArray("CATTAATTCTTT"), new Sequence_ByteArray("A-TTAATT-TTT"));
            alignment.AlignmentFormater = new PairwiseAlignmentFormater_WaterEmboss();
            Assert.AreEqual(alignment[0], new char[] { 'C', '.', 'A' });
            Assert.AreEqual(alignment[1], new char[] { 'A', ' ', '-' });
        }

        [Test]
        public void Test_PairwiseAlignmentBuilder()
        {
            PairwiseAlignmentBuilder ab = new PairwiseAlignmentBuilder();
            ab.Append_3_prime('A', 'A');
            ab.Append_3_prime('T', 'T');
            ab.Append_5_prime('C', 'C');
            ab.Append_5_prime('G', 'G');
            PairwiseAlignment al = ab.GetAlignment();

            Assert.AreEqual(al.DatabaseSequence.ToString(), "GCAT");
            Assert.AreEqual(al.QuerySequence.ToString(), "GCAT");
            Assert.AreEqual(al.SimilaritySequence.ToString(), "||||");
            ab.Append_3_prime('T', 'A');
            ab.Append_5_prime('C', 'G');
            al = ab.GetAlignment();
            Assert.AreEqual(al.DatabaseSequence.ToString(), "CGCATT");
            Assert.AreEqual(al.QuerySequence.ToString(), "GGCATA");
            Assert.AreEqual(al.SimilaritySequence.ToString(), ".||||.");
            al = ab.GetAlignment();
            Assert.AreEqual(al.DatabaseSequence.ToString(), "CGCATT");
            Assert.AreEqual(al.QuerySequence.ToString(), "GGCATA");
            Assert.AreEqual(al.SimilaritySequence.ToString(), ".||||.");

            ///Rount Two
            ab.Append_3_prime(new PairwiseAlignment(new Sequence_ByteArray("ACT"), new Sequence_ByteArray("ACT")));
            ab.Append_5_prime(new PairwiseAlignment(new Sequence_ByteArray("GTC"), new Sequence_ByteArray("GTC")));

            al = ab.GetAlignment();
            Assert.AreEqual(al.DatabaseSequence.ToString(), "GTCCGCATTACT");
            Assert.AreEqual(al.QuerySequence.ToString(), "GTCGGCATAACT");
            Assert.AreEqual(al.SimilaritySequence.ToString(), "|||.||||.|||");

            ab = new PairwiseAlignmentBuilder(al);
            ab.Append_3_prime('C', 'C');
            ab.Append_5_prime('A', 'A');
            al = ab.GetAlignment();
            Assert.AreEqual(al.DatabaseSequence.ToString(), "AGTCCGCATTACTC");
            Assert.AreEqual(al.QuerySequence.ToString(), "AGTCGGCATAACTC");
            Assert.AreEqual(al.SimilaritySequence.ToString(), "||||.||||.||||");

        }

        [Test]
        public void Test_SmithWatermanGotoh()
        {
            ISequenceContainer database = new Sequence_ByteArray("TTTTTTTTTTACGTACGTACGTTTTTTTTTTT");
            ISequenceContainer query = new Sequence_ByteArray("CCCCCCCCCCACGTCCCCACGTCCCCCCCCCC");
            SmithWatermanGotoh sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGT\r\n||||.|..||||\r\nACGTCCCCACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Score, 6);
            Assert.AreEqual(sml.Start_Query, 11);
            Assert.AreEqual(sml.End_Database, 22);
            Assert.AreEqual(sml.End_Query, 22);

            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1, 0.1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGT----ACGT\r\n||||    ||||\r\nACGTCCCCACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Start_Query, 11);
            Assert.AreEqual(sml.End_Database, 18);
            Assert.AreEqual(sml.End_Query, 22);
            Assert.IsTrue(Math.Abs(sml.Score - 6.7F) < 0.01);


            database = new Sequence_ByteArray("TTTTTT");
            query = new Sequence_ByteArray("AAAAAA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.Score, 0);

            database = new Sequence_String("TTTTTT");
            query = new Sequence_String("TTTTTT");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTT\r\n||||||\r\nTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 6);
            Assert.AreEqual(sml.End_Query, 6);
            Assert.AreEqual(sml.Score, 6.0);

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("A");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 1);


            database = new Sequence_ByteArray("TTTTA");
            query = new Sequence_ByteArray("CCCCA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 5);
            Assert.AreEqual(sml.Start_Query, 5);
            Assert.AreEqual(sml.End_Database, 5);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 1);


            database = new Sequence_ByteArray("ATTTT");
            query = new Sequence_ByteArray("ACCCT");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 1);

            database = new Sequence_ByteArray("ATTTT");
            query = new Sequence_ByteArray("ACCTT");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TT\r\n||\r\nTT\r\n");
            Assert.AreEqual(sml.Start_Database, 2);
            Assert.AreEqual(sml.Start_Query, 4);
            Assert.AreEqual(sml.End_Database, 3);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 2);

            database = new Sequence_ByteArray("AATTCCGGAATTTTCCGGAATTTCGG");
            query = new Sequence_ByteArray("AATTCCCCGGAATTCCGGAATTCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AATT--CCGGAATTTTCCGGAATTTCGG\r\n||||  ||||||  ||||||||||.|||\r\nAATTCCCCGGAA--TTCCGGAATTCCGG\r\n");
            Assert.AreEqual(sml.Score, 16);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 26);
            Assert.AreEqual(sml.End_Query, 26);

            database = new Sequence_String("AATTCCGGAATTTTCCGGAATTTCGG");
            query = new Sequence_String("AATTCCCCGGAATTCCGGAATTCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AATT--CCGGAATTTTCCGGAATTTCGG\r\n||||  ||||||  ||||||||||.|||\r\nAATTCCCCGGAA--TTCCGGAATTCCGG\r\n");
            Assert.AreEqual(sml.Score, 16);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 26);
            Assert.AreEqual(sml.End_Query, 26);

            //Score test match
            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("A");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2F, 1F, 2, 1F));
            Assert.AreEqual(sml.Score, 2);

            database = new Sequence_ByteArray("ATA");
            query = new Sequence_ByteArray("AAA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2F, 1F, 2, 1F));
            Assert.AreEqual(sml.Score, 3);

            database = new Sequence_ByteArray("AATTAATTAATTACCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 21);

            database = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTACCGGCCGGCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 21);

            database = new Sequence_ByteArray("AATTAATTAATTAACCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 20);

            database = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTAACCGGCCGGCCGG");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 20);

            database = new Sequence_ByteArray("ANA");
            query = new Sequence_ByteArray("AAA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 2);
            Assert.AreEqual(al.ToString(), "ANA\r\n| |\r\nAAA\r\n");

            database = new Sequence_ByteArray("TTTNXAAA");
            query = new Sequence_ByteArray("TTTAAAAA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 5);
            Assert.AreEqual(al.ToString(), "TTTNXAAA\r\n|||  |||\r\nTTTAAAAA\r\n");

            database = new Sequence_ByteArray("ANA");
            query = new Sequence_ByteArray("ANA");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 2);
            Assert.AreEqual(al.ToString(), "ANA\r\n| |\r\nANA\r\n");


            database = new Sequence_ByteArray("CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTAATAGCAAAGTTCTTCTACTTCGCAAA");
            query = new Sequence_ByteArray("GTATCAGATTTAAAGCTCCATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTA-----ATAGCAAAGTTCTTCTACTTCGC\r\n||||||.||||||||||    ||.||||||||.|||||||||||||||||.||.||.||.||||||||||.||||||     .|.||||||||||||||||||||\r\nTATCAGATTTAAAGCTC----CATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGC\r\n");

            Assert.AreEqual(sml.Score, 59);
            Assert.AreEqual(sml.Start_Database, 8);
            Assert.AreEqual(sml.Start_Query, 2);
            Assert.AreEqual(sml.End_Database, 107);
            Assert.AreEqual(sml.End_Query, 102);


            database = new Sequence_ByteArray("CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTAATAGCAAAGTTCTTCTACTTCGCAAA");
            query = new Sequence_ByteArray("GTATCAGATTTAAAGCTCCATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT");
            sml = new SmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion());
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTA-----ATAGCAAAGTTCTTCTACTTCGC\r\n||||||.||||||||||    ||.||||||||.|||||||||||||||||.||.||.||.||||||||||.||||||     .|.||||||||||||||||||||\r\nTATCAGATTTAAAGCTC----CATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGC\r\n");

            Assert.AreEqual(sml.Score, 97.56012);
            Assert.AreEqual(sml.Start_Database, 8);
            Assert.AreEqual(sml.Start_Query, 2);
            Assert.AreEqual(sml.End_Database, 107);
            Assert.AreEqual(sml.End_Query, 102);


#if (DEBUG && FILE) 

            StreamWriter sw = new StreamWriter("G:\\Temp\\test.txt", false, Encoding.ASCII);
            sw.Write(matrix);
            sw.Close();
#endif
        }

        [Test]
        public void Test_NeedlemanWunschGotoh()
        {


            ISequenceContainer database = new Sequence_ByteArray("TTTTTTTTTTACGTACGTACGTTTTTTTTTTT");
            ISequenceContainer query = new Sequence_ByteArray("CCCCCCCCCCACGTCCCCACGTCCCCCCCCCC");
            NeedlemanWunschGotoh sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTTACGTACGTACGTTTTTTTTTTT\r\n..........||||.|..||||..........\r\nCCCCCCCCCCACGTCCCCACGTCCCCCCCCCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, -14.0);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 32);
            Assert.AreEqual(sml.End_Query, 32);


            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1, 0.1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "----------TTTTTTTTTTACGT----ACGTA--------CGTTTTTTTTTTT\r\n                    ||||    ||||.        |            \r\nCCCCCCCCCC----------ACGTCCCCACGTCCCCCCCCCC------------\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 32);
            Assert.AreEqual(sml.End_Query, 32);
            Assert.IsTrue(Math.Abs(sml.Score + 0.9F) < 0.01);

            sml = new NeedlemanWunschGotoh(new Sequence_ByteArray("A"), new Sequence_ByteArray("T"), SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1, 0.1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n.\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, -1);

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("A");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1, 0.1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 1);

            database = new Sequence_ByteArray("CCCCAAA");
            query = new Sequence_ByteArray("AAA");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCCCAAA\r\n    |||\r\n----AAA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, -5);

            database = new Sequence_ByteArray("AAACCCC");
            query = new Sequence_ByteArray("AAA");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAACCCC\r\n|||    \r\nAAA----\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, -5);

            database = new Sequence_ByteArray("ATACCCC");
            query = new Sequence_ByteArray("AAA");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ATACCCC\r\n|.|    \r\nAAA----\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, -7);

            database = new Sequence_ByteArray("CCCCAAACCCC");
            query = new Sequence_ByteArray("AAA");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCCCAAACCCC\r\n    |||    \r\n----AAA----\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 11);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, -13);

            database = new Sequence_ByteArray("AAA");
            query = new Sequence_ByteArray("CCCCAAACCCC");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "----AAA----\r\n    |||    \r\nCCCCAAACCCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 3);
            Assert.AreEqual(sml.End_Query, 11);
            Assert.AreEqual(sml.Score, -13);


            database = new Sequence_ByteArray("CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTAATAGCAAAGTTCTTCTACTTCGCAAA");
            query = new Sequence_ByteArray("GTATCAGATTTAAAGCTCCATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTA-----ATAGCAAAGTTCTTCTACTTCGC--------AAA\r\n      .||||||.||||||||||    ||.||||||||.|||||||||||||||||.||.||.||.||||||||||.||||||     .|.||||||||||||||||||||        ...\r\n------GTATCAGATTTAAAGCTC----CATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Score, 33.0F);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 110);
            Assert.AreEqual(sml.End_Query, 113);


            database = new Sequence_ByteArray("CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTAATAGCAAAGTTCTTCTACTTCGCAAA");
            query = new Sequence_ByteArray("GTATCAGATTTAAAGCTCCATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT");
            sml = new NeedlemanWunschGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion());
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTA-----ATAGCAAAGTTCTTCTACTTCGC--------AAA\r\n      .||||||.|||||||||||.|    ||||||||.|||||||||||||||||.||.||.||.||||||||||.||||||     .|.||||||||||||||||||||        ...\r\n------GTATCAGATTTAAAGCTCCAT----GGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Score, 53.7200317f);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 110);
            Assert.AreEqual(sml.End_Query, 113);



        }

        [Test]
        public void Test_SmithWatermanGotoh_454P()
        {


            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            SmithWatermanGotoh_454P sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 2));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "99");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "99");
            Assert.AreEqual(sml.ToString(), "\t$\tT\t\n$\t0,00\t0,00\nT\t0,00\t2,00\n");

            database = new Sequence_ByteArray("TTT");
            query = new Sequence_ByteArray("TTTATATAT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.End_Database, 3);
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTTTTTTF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "9951");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "9951999999");

            database = new Sequence_ByteArray("TTTTTTTTTT");
            query = new Sequence_ByteArray("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 10);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FFFFFFFFFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFFFFFFFFFF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "99876544321");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "998766543321");

            database = new Sequence_String("TTTTTTTTTT");
            query = new Sequence_String("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 10);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FFFFFFFFFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFFFFFFFFFF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "99876544321");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "998766543321");

            database = new Sequence_ByteArray("CGTTTT");
            query = new Sequence_ByteArray("AAATTT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FTTFFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTFFF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "9999765");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "9975975");

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("T");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);


            ///
            ///Test for roberts law
            ///

            database = new Sequence_ByteArray("CGATTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTT-CGTCGT\r\n|||||| ||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTCGTCGT\r\n||||||.||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");



            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 20, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");


            ///
            ///
            ///Unspecific
            ///


            database = new Sequence_ByteArray("ACGTACGTTTTTTTACGTACGTACGAAAAAAATCGTCGTGTA");
            query = new Sequence_ByteArray("ACGAACGTTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTTTTTTTACGTACGTACGAAAAAAA-TCGTCGTGTA\r\n|||.||||||||| |||||||.|||||||||| ||||||.|||\r\nACGAACGTTTTTT-ACGTACGCACGAAAAAAAATCGTCGAGTA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 42);
            Assert.AreEqual(sml.End_Query, 42);
            Assert.IsTrue(Math.Abs(sml.Score - 47.14) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CCCCCC-GGGGGGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| |||||| ||||||| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 41);
            Assert.AreEqual(sml.End_Query, 47);
            Assert.IsTrue(Math.Abs(sml.Score - 106.0F) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CTCCCCCCAGAGAGG--AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||.|||||.|.|.||  ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");


            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_String("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_String("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");



            database = new Sequence_ByteArray("CGTCGTTATAAATTAAAA");
            query = new Sequence_ByteArray("CGTCGTTAAATAAAAAAA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTTA--TAAATTAAAA\r\n||||||||  ||||  ||||\r\nCGTCGTTAAATAAA--AAAA\r\n");



            database = new Sequence_ByteArray("CGTCGATTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGAATTTCGTCGT");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGA-TTTCGTCGT\r\n|||||| |||||||||\r\nCGTCGAATTTCGTCGT\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");



        }

        [Test]
        public void Test_SmithWatermanGotoh_454M()
        {


            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            SmithWatermanGotoh_454M sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 2));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "99");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "99");
            //Assert.AreEqual(sml.ToString(), "T\t$\t\t\n$\t2,00\t0,00\nT\t0,00\t0,00\n");

            database = new Sequence_ByteArray("TTT");
            query = new Sequence_ByteArray("TTTATATAT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTTTTTTF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "1599");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "1599999999");

            database = new Sequence_String("TTTTTT");
            query = new Sequence_String("TTTTTT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTT\r\n||||||\r\nTTTTTT\r\n");


            database = new Sequence_ByteArray("CGTTTT");
            query = new Sequence_ByteArray("AAATTT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FTTFFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTFFF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "9914699");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "1591599");

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("T");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);


            ///
            ///Test for roberts law
            ///

            database = new Sequence_ByteArray("CGATTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTT-CGTCGT\r\n|||||| ||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTCGTCGT\r\n||||||.||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n||| |||.||||||\r\nCGA-TTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n||| |||.||||||\r\nCGA-TTTACGTCGT\r\n");



            database = new Sequence_ByteArray("ACACAAAAACAAAACAAAGACAACAAAAC");
            query = new Sequence_ByteArray("ACACAAAAACAAAACAAATGACAACAAAAC");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACACAAAAACAAAACAAA-GACAACAAAAC\r\n|||||||||||||||||| |||||||||||\r\nACACAAAAACAAAACAAATGACAACAAAAC\r\n");


            database = new Sequence_ByteArray("ACACAAAAACAAAACAAAGACAACAAAAC");
            query = new Sequence_ByteArray("ACACAAAAACAAAACAAATGACAACAAAAC");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACACAAAAACAAAACAAA-GACAACAAAAC\r\n|||||||||||||||||| |||||||||||\r\nACACAAAAACAAAACAAATGACAACAAAAC\r\n");



            database = new Sequence_ByteArray("ACACAAAAACAAAACAAAGACAACAAAAC");
            query = new Sequence_ByteArray("ACACAAAAACAAAACAAATGACAACAAAAC");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 20, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACACAAAAACAAAACAAA-GACAACAAAAC\r\n|||||||||||||||||| |||||||||||\r\nACACAAAAACAAAACAAATGACAACAAAAC\r\n");


            ///
            ///
            ///Unspecific
            ///


            database = new Sequence_ByteArray("ACGTACGTTTTTTTACGTACGTACGAAAAAAATCGTCGTGTA");
            query = new Sequence_ByteArray("ACGAACGTTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTTTTTTTACGTACGTACG-AAAAAAATCGTCGTGTA\r\n|||.||| |||||||||||||.||| |||||||||||||.|||\r\nACGAACG-TTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 42);
            Assert.AreEqual(sml.End_Query, 42);
            Assert.IsTrue(Math.Abs(sml.Score - 47.14) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGT--AAAA-TTTTT-CCCCCC-GGGGGGG-AAAAAAACGTCGT\r\n||||||  |||| ||||| |||||| ||||||| |||||||||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 41);
            Assert.AreEqual(sml.End_Query, 47);
            Assert.IsTrue(Math.Abs(sml.Score - 106.00F) < 0.01);

            database = new Sequence_String("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_String("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGT--AAAA-TTTTT-CCCCCC-GGGGGGG-AAAAAAACGTCGT\r\n||||||  |||| ||||| |||||| ||||||| |||||||||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 41);
            Assert.AreEqual(sml.End_Query, 47);
            Assert.IsTrue(Math.Abs(sml.Score - 106.00F) < 0.01);





            database = new Sequence_ByteArray("CGTCGTTATAAATTAAAA");
            query = new Sequence_ByteArray("CGTCGTTAAATAAAAAAA");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTT--ATAAATTAAAA\r\n|||||||  |||||  ||||\r\nCGTCGTTAAATAAA--AAAA\r\n");






            database = new Sequence_ByteArray("TAATTTAAATATAAAATTTCACTTAATATA");
            query = new Sequence_ByteArray("TAATTTAAATATTAAAATTTCACTTAATATA");
            sml = new SmithWatermanGotoh_454M(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TAATTTAAATA-TAAAATTTCACTTAATATA\r\n||||||||||| |||||||||||||||||||\r\nTAATTTAAATATTAAAATTTCACTTAATATA\r\n");


        }

        [Test]
        public void Test_PairwiseNucleotideSequenceAlignment()
        {
            PairwiseAlignment a1 = new PairwiseAlignment(new Sequence_ByteArray("AAAAACCCCCTTTTT-GGGG"), new Sequence_ByteArray("AAAA-CCCC-TTTGTGGGGG"));
            PairwiseNucleotideSequenceAlignment pw1 = new PairwiseNucleotideSequenceAlignment(a1, "ferdinand", "query1", 200, 10, 300, 110);
            pw1.PlusPlusStrand = true;
            pw1.Significator = new BlastSignificator_Score(10);
            pw1.Score = 10;

            Assert.AreEqual(pw1.DatabaseParent, "ferdinand");
            Assert.AreEqual(pw1.QueryParent, "query1");
            Assert.AreEqual(pw1.StartDatabase, 200);
            Assert.AreEqual(pw1.Start, 200);
            Assert.AreEqual(pw1.EndDatabase, 300);
            Assert.AreEqual(pw1.End, 300);
            Assert.AreEqual(pw1.EndQuery, 110);
            Assert.AreEqual(pw1.LengthAlignmentWithGaps, 20);
            Assert.AreEqual(pw1.Length, 101); ///Length of the Database
            Assert.AreEqual(pw1.LengthAlignmentWithoutGaps, 17);
            Assert.AreEqual(pw1.LengthAlignedDatabase, 101);
            Assert.AreEqual(pw1.LengthAlignedQuery, 101);
            Assert.AreEqual(pw1.PlusPlusStrand, true);
            Assert.AreEqual(pw1.Significator.GetType(), (new BlastSignificator_Score(7)).GetType());
            Assert.AreEqual(pw1.IsSignificant(), true);
            Assert.AreEqual(pw1.Gaps, 3);
            Assert.AreEqual(pw1.Hits, 16);
            Assert.AreEqual(pw1.IsRoot, false);
            Assert.AreEqual(pw1.Score, 10);
            Assert.AreEqual(pw1.Count_LongGaps_Database(2), 0);
            Assert.AreEqual(pw1.Count_LongGaps_Database(1), 1);
            Assert.AreEqual(pw1.Count_LongGaps_Query(2), 0);
            Assert.AreEqual(pw1.Count_LongGaps_Query(1), 2);
            Assert.AreEqual(pw1.Count_LongGaps(2), 0);
            Assert.AreEqual(pw1.Count_LongGaps(1), 3);
            Assert.AreEqual(pw1.Transitions, 0);
            Assert.AreEqual(pw1.Transversions, 1);


            a1 = new PairwiseAlignment(new Sequence_ByteArray("AAAAAAGGGGGGGCCC---AAA"), new Sequence_ByteArray("AAA---AAATTTTCCCAAAAAA"));
            PairwiseNucleotideSequenceAlignment pw2 = new PairwiseNucleotideSequenceAlignment(a1, "ferdinand", "query1", 400, 210, 500, 310);
            pw2.PlusPlusStrand = true;
            Assert.AreEqual(pw2.Count_LongGaps(3), 2);
            Assert.AreEqual(pw2.Count_LongGaps(4), 0);
            Assert.AreEqual(pw2.Transversions, 4);
            Assert.AreEqual(pw2.Transitions, 3);
            Assert.AreEqual(pw2.DatabaseParent, "ferdinand");
            Assert.AreEqual(pw2.QueryParent, "query1");

            List<PairwiseNucleotideSequenceAlignment> lp = new List<PairwiseNucleotideSequenceAlignment>();
            lp.Add(pw1);
            lp.Add(pw2);
            CompositePairwiseNucleotideSequenceAlignment cpa = new CompositePairwiseNucleotideSequenceAlignment(lp);
            Assert.AreEqual(cpa.DatabaseParent, "ferdinand");
            Assert.AreEqual(cpa.QueryParent, "query1");
            Assert.AreEqual(cpa.StartDatabase, 200);
            Assert.AreEqual(cpa.StartQuery, 10);
            Assert.AreEqual(cpa.EndDatabase, 500);
            Assert.AreEqual(cpa.EndQuery, 310);
            Assert.AreEqual(cpa.Count_LongGaps(99), 2);
            Assert.AreEqual(cpa.Count_LongGaps(100), 0);
            Assert.AreEqual(cpa.Count_LongGaps_Database(1), 3);
            Assert.AreEqual(cpa.Count_LongGaps_Query(1), 4);
            Assert.AreEqual(cpa.Transitions, 3);
            Assert.AreEqual(cpa.Transversions, 5);
            Assert.AreEqual(cpa.Count_SubAlignments, 2);
            Assert.AreEqual(cpa[0], pw1);
            Assert.AreEqual(cpa[1], pw2);
            Assert.AreEqual(cpa.Gaps, 9);


            pw1 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AATA-"), new Sequence_ByteArray("AAAAA")), "ferdinand", "query1", 200, 10, 300, 110);

            Assert.AreEqual(pw1.SimilarityWithGaps, 60.0F);
            Assert.AreEqual(pw1.SimilarityWithoutGaps, 75.0F);

            ///
            ///TEST Subalignment
            ///
            a1 = new PairwiseAlignment(new Sequence_ByteArray("AAA---AAAGGGGGGGCCC---AAA"), new Sequence_ByteArray("AAAAAA---AAATTTTCCCAAAAAA"));
            PairwiseNucleotideSequenceAlignment pw3 = new PairwiseNucleotideSequenceAlignment(a1, "ferdinand", "query1", 400, 210, 500, 310);
            pw3.PlusPlusStrand = true;
            PairwiseNucleotideSequenceAlignment pw4 = pw3.SubAlignment_RelativeToQuery(212);
            Assert.AreEqual(pw4.StartQuery, 212);
            Assert.AreEqual(pw4.StartDatabase, 402);
            Assert.AreEqual(pw4.EndQuery, 310);
            Assert.AreEqual(pw4.EndDatabase, 500);
            Assert.AreEqual(pw4.Alignment.DatabaseSequence.ToString(), "A---AAAGGGGGGGCCC---AAA");
            Assert.AreEqual(pw4.Alignment.QuerySequence.ToString(), "AAAA---AAATTTTCCCAAAAAA");
            Assert.AreEqual(pw4.PlusPlusStrand, true);

            pw4 = pw3.SubAlignment_RelativeToQuery(213);
            Assert.AreEqual(pw4.StartQuery, 216);
            Assert.AreEqual(pw4.StartDatabase, 406);
            Assert.AreEqual(pw4.EndQuery, 310);
            Assert.AreEqual(pw4.EndDatabase, 500);
            Assert.AreEqual(pw4.Alignment.DatabaseSequence.ToString(), "GGGGGGGCCC---AAA");
            Assert.AreEqual(pw4.Alignment.QuerySequence.ToString(), "AAATTTTCCCAAAAAA");
            Assert.AreEqual(pw4.PlusPlusStrand, true);

        }
    }

    [TestFixture]
    public class Aln_Misc
    {

        public static void Start()
        {
            Aln_Misc p = new Aln_Misc();
            p.Test_SmithWaterman_454P();
            p.Test_SmithWaterman();
            p.Test_SmithWatermanGotoh_DynamicBanded_3P();
            p.Test_SmithWatermanGotoh_DynamicBanded_5P();
            p.Test_SmitWatermanGotoh_DynamicBanded_454_3p();
            p.Test_SmitWatermanGotoh_DynamicBanded_454_5p();
        }

        [Test]
        public void Test_SmithWaterman()
        {
            ///Implementation of the original smith waterman algorithm, should be extremly time consuming but it is a highly reliable implementation

            ISequenceContainer database = new Sequence_ByteArray("TTTTTTTTTTACGTACGTACGTTTTTTTTTTT");
            ISequenceContainer query = new Sequence_ByteArray("CCCCCCCCCCACGTCCCCACGTCCCCCCCCCC");
            SmithWaterman sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5F, 1F));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGT\r\n||||.|..||||\r\nACGTCCCCACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Score, 6);
            Assert.AreEqual(sml.Start_Query, 11);
            Assert.AreEqual(sml.End_Database, 22);
            Assert.AreEqual(sml.End_Query, 22);

            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1, 0.1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGT----ACGT\r\n||||    ||||\r\nACGTCCCCACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Start_Query, 11);
            Assert.AreEqual(sml.End_Database, 18);
            Assert.AreEqual(sml.End_Query, 22);
            Assert.AreEqual(sml.Score, 6.7F);


            database = new Sequence_ByteArray("TTTTTT");
            query = new Sequence_ByteArray("AAAAAA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.Score, 0);

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("A");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 1);


            database = new Sequence_ByteArray("TTTTA");
            query = new Sequence_ByteArray("CCCCA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 5);
            Assert.AreEqual(sml.Start_Query, 5);
            Assert.AreEqual(sml.End_Database, 5);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 1);


            database = new Sequence_ByteArray("ATTTT");
            query = new Sequence_ByteArray("ACCCT");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 1);

            database = new Sequence_ByteArray("ATTTT");
            query = new Sequence_ByteArray("ACCTT");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TT\r\n||\r\nTT\r\n");
            Assert.AreEqual(sml.Start_Database, 2);
            Assert.AreEqual(sml.Start_Query, 4);
            Assert.AreEqual(sml.End_Database, 3);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 2);

            database = new Sequence_ByteArray("AATTCCGGAATTTTCCGGAATTTCGG");
            query = new Sequence_ByteArray("AATTCCCCGGAATTCCGGAATTCCGG");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AATT--CCGGAATTTTCCGGAATTTCGG\r\n||||  ||||||  ||||||||||.|||\r\nAATTCCCCGGAA--TTCCGGAATTCCGG\r\n");
            Assert.AreEqual(sml.Score, 16);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 26);
            Assert.AreEqual(sml.End_Query, 26);

            //Score test match
            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("A");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2F, 1F, 2, 1F));
            Assert.AreEqual(sml.Score, 2);

            database = new Sequence_ByteArray("ATA");
            query = new Sequence_ByteArray("AAA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2F, 1F, 2, 1F));
            Assert.AreEqual(sml.Score, 3);

            database = new Sequence_ByteArray("AATTAATTAATTACCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 21);

            database = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTACCGGCCGGCCGG");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 21);

            database = new Sequence_ByteArray("AATTAATTAATTAACCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 20);

            database = new Sequence_ByteArray("AATTAATTAATTCCGGCCGGCCGG");
            query = new Sequence_ByteArray("AATTAATTAATTAACCGGCCGGCCGG");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 20);

            database = new Sequence_ByteArray("ANA");
            query = new Sequence_ByteArray("AAA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 2);
            Assert.AreEqual(al.ToString(), "ANA\r\n| |\r\nAAA\r\n");

            database = new Sequence_ByteArray("TTTNXAAA");
            query = new Sequence_ByteArray("TTTAAAAA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 5);
            Assert.AreEqual(al.ToString(), "TTTNXAAA\r\n|||  |||\r\nTTTAAAAA\r\n");

            database = new Sequence_ByteArray("ANA");
            query = new Sequence_ByteArray("ANA");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 3, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.Score, 2);
            Assert.AreEqual(al.ToString(), "ANA\r\n| |\r\nANA\r\n");


            database = new Sequence_ByteArray("CCCCCCCTATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTAATAGCAAAGTTCTTCTACTTCGCAAA");
            query = new Sequence_ByteArray("GTATCAGATTTAAAGCTCCATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGCTTTTTTTTTTT");
            sml = new SmithWaterman(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1F));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATCAGGTTTAAAGCTCCGTGCACGGAAGATATTTAATTCTCAACGCGACTTGCGTTGAGAATTTTCCACAGAGTTA-----ATAGCAAAGTTCTTCTACTTCGC\r\n||||||.||||||||||    ||.||||||||.|||||||||||||||||.||.||.||.||||||||||.||||||     .|.||||||||||||||||||||\r\nTATCAGATTTAAAGCTC----CATGGAAGATAGTTAATTCTCAACGCGACATGTGTAGAAAATTTTCCACGGAGTTAGACGTGTTGCAAAGTTCTTCTACTTCGC\r\n");
            Assert.AreEqual(sml.Score, 59);
            Assert.AreEqual(sml.Start_Database, 8);
            Assert.AreEqual(sml.Start_Query, 2);
            Assert.AreEqual(sml.End_Database, 107);
            Assert.AreEqual(sml.End_Query, 102);







        }

        [Test]
        public void Test_SmithWaterman_454P()
        {



            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            SmithWaterman_454P sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 2));
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "99");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "99");
            Assert.AreEqual(sml.ToString(), "\t$\tT\t\n$\t0,00\t0,00\nT\t0,00\t2,00\n");

            database = new Sequence_ByteArray("TTT");
            query = new Sequence_ByteArray("TTTATATAT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTTTTTTF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "9976");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "9976999999");

            database = new Sequence_ByteArray("CGTTTT");
            query = new Sequence_ByteArray("AAATTT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FTTFFFF");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FFFTFFF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "9999764");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "9976976");

            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("T");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);


            database = new Sequence_ByteArray("TTTTTTTTTTACGTACGTACGTTTTTTTTTTT");
            query = new Sequence_ByteArray("CCCCCCCCCCACGTCCCCACGTCCCCCCCCCC");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGT\r\n||||.|..||||\r\nACGTCCCCACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Score, 6);
            Assert.AreEqual(sml.Start_Query, 11);
            Assert.AreEqual(sml.End_Database, 22);
            Assert.AreEqual(sml.End_Query, 22);
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "554433222115555555555554433322211");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "554433222115555543355555443322211");


            database = new Sequence_ByteArray("ACGTACGTTTTTTTACGTACGTACGAAAAAAATCGTCGTGTA");
            query = new Sequence_ByteArray("ACGAACGTTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTTTTTTTACGTACGTACGAAAAAAA-TCGTCGTGTA\r\n|||.||||||||| |||||||.|||||||||| ||||||.|||\r\nACGAACGTTTTTT-ACGTACGCACGAAAAAAAATCGTCGAGTA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 42);
            Assert.AreEqual(sml.End_Query, 42);
            Assert.IsTrue(Math.Abs(sml.Score - 52.95F) < 0.01);

            database = new Sequence_ByteArray("ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT");
            query = new Sequence_ByteArray("ACGTACGTACGTTTTTGTTTTACGTACGTACGTACGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT\r\n||||||||||||||||.||||  ||||||||||||||||\r\nACGTACGTACGTTTTTGTTTT--ACGTACGTACGTACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 39);
            Assert.AreEqual(sml.End_Query, 37);
            Assert.IsTrue(Math.Abs(sml.Score - 54.63F) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CCCCCC-GGGGGGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| |||||| ||||||| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 41);
            Assert.AreEqual(sml.End_Query, 47);
            Assert.IsTrue(Math.Abs(sml.Score - 56.25F) < 0.01);
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "999999997649764397643297654329765432999999");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "999999997643297643297654329876543298765432999999");

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            al.AlignmentFormater = new PairwiseAlignmentFormater_Blast();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CTCCCC---AGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| || |||    | | || ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTTTCTTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGTTTTCTTTTCGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTTT-CTTT-CGTCGT\r\n|||||||| |||| ||||||\r\nCGTCGTTTTCTTTTCGTCGT\r\n");


            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            al.AlignmentFormater = new PairwiseAlignmentFormater_Blast();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CTCCCC---AGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| || |||    | | || ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            al.AlignmentFormater = new PairwiseAlignmentFormater_Blast();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n||||||||||  ||||| || |||    | | || ||||||| ||||||\r\nCGTCGTAAAA--TTTTT-CTCCCC---AGAGAGG-AAAAAAA-CGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTTTTTTTTTTTTTTTTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGTATATATATATTTTTCGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            al.AlignmentFormater = new PairwiseAlignmentFormater_Blast();
            Assert.AreEqual(al.ToString(), "CGTCGTTTTTTTTTTTTTTTTTCGTCGT\r\n|||||| | | | | |||||  ||||||\r\nCGTCGTATATATATATTTTT--CGTCGT\r\n");


            database = new Sequence_ByteArray("CGTCGTCGTTTTATTTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGTCGTATATATATATTTTTCGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 3, 9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(sml.ToStringVariableGapPenaltyDatabase(), "999999999976499764999999");
            Assert.AreEqual(sml.ToStringBoundariesDatabase(), "FTTTTTTTTFFFTTFFFTTTTTTF");
            Assert.AreEqual(sml.ToStringVariableGapPenaltyQuery(), "999999999999999999997643999999");
            Assert.AreEqual(sml.ToStringBoundariesQuery(), "FTTTTTTTTTTTTTTTTTTFFFFTTTTTTF");

            Assert.AreEqual(al.ToString(), "CGTCGTCGT-----TTTATTTT-CGTCGT\r\n|||||||||     |.|||||| ||||||\r\nCGTCGTCGTATATATATATTTTTCGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTTTCTTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGTTTTCTTTTCGTCGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTTT-CTTT-CGTCGT\r\n|||||||| |||| ||||||\r\nCGTCGTTTTCTTTTCGTCGT\r\n");



            database = new Sequence_ByteArray("ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT");
            query = new Sequence_ByteArray("ACGTACGTACGTTTTGTTTTTACGTACGTACGTACGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();

            Assert.AreEqual(al.ToString(), "ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT\r\n|||||||||||||||.|||||  ||||||||||||||||\r\nACGTACGTACGTTTTGTTTTT--ACGTACGTACGTACGT\r\n");



            database = new Sequence_ByteArray("ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT");
            query = new Sequence_ByteArray("ACGTACGTACGTTTTTTGTTTTTACGTACGTACGTACGT");
            sml = new SmithWaterman_454P(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1));
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGTTTTTTTTTTTTACGTACGTACGTACGT\r\n|||||||||||||||||.|||||||||||||||||||||\r\nACGTACGTACGTTTTTTGTTTTTACGTACGTACGTACGT\r\n");






        }

        [Test]
        public void Test_SmithWatermanGotoh_DynamicBanded_3P()
        {
            //
            // 3'-Prime version
            //
            ISequenceContainer database = new Sequence_ByteArray("ATTCCTTAA");
            ISequenceContainer query = new Sequence_ByteArray("TTTCCTCC");
            SmithWatermanGotoh_DynamicBanded_3p sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5F, 1.0F), 0);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ATTCCT\r\n.|||||\r\nTTTCCT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 4);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 6);
            Assert.AreEqual(sml.End_Query, 6);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("AAAATTTTT");
            query = new Sequence_ByteArray("TTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAAATTTTT\r\n    |||||\r\n----TTTTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 9);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 3.7F);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("AAAATTTTT");
            query = new Sequence_ByteArray("TTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 2);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAAATTT\r\n    |||\r\n----TTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, 1.7F);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_String("TTTTT");
            query = new Sequence_String("TTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F),0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTT\r\n|||||\r\nTTTTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 5);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 5);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("C");
            query = new Sequence_ByteArray("G");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 2);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.EndofDynamicExtension, true);
            Assert.AreEqual(sml.Score, 0.0F);


            database = new Sequence_ByteArray("AAAATTTTT");
            query = new Sequence_ByteArray("TTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 2);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAAATTT\r\n    |||\r\n----TTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, 1.7F);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("TTTCCAAAATTTTT");
            query = new Sequence_ByteArray("TTTCCTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 5);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTCC\r\n|||||\r\nTTTCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 5);
            Assert.AreEqual(sml.End_Query, 5);
            Assert.AreEqual(sml.Score, 5.0F);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("TTTCCAAAATTTTT");
            query = new Sequence_ByteArray("TTTCCTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTCCAAAAT\r\n|||||    |\r\nTTTCC----T\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.End_Query, 6);
            Assert.IsTrue(Math.Abs(sml.Score - 4.7F) < 0.0001);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_String("TTTCCAAAATTTTT");
            query = new Sequence_String("TTTCCTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTCCAAAAT\r\n|||||    |\r\nTTTCC----T\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.End_Query, 6);
            Assert.IsTrue(Math.Abs(sml.Score - 4.7F) < 0.0001);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


        }

        [Test]
        public void Test_SmithWatermanGotoh_DynamicBanded_5P()
        {
            ISequenceContainer database;
            ISequenceContainer query;
            PairwiseAlignment al;

            database = new Sequence_ByteArray("AAAATTTTT");
            query = new Sequence_ByteArray("TTTTT");
            SmithWatermanGotoh_DynamicBanded_5p sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 0);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTT\r\n|||||\r\nTTTTT\r\n");
            Assert.AreEqual(sml2.Start_Database, -5);
            Assert.AreEqual(sml2.Start_Query, -5);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.AreEqual(sml2.Score, 5.0F);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("TTTTTAAAA");
            query = new Sequence_ByteArray("TTTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 0);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTAAAA\r\n|||||    \r\nTTTTT----\r\n");
            Assert.AreEqual(sml2.Start_Database, -9);
            Assert.AreEqual(sml2.Start_Query, -5);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml2.Score - 3.7F) < 0.001);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);

            database = new Sequence_String("TTTTT");
            query = new Sequence_String("TTTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 0);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTT\r\n|||||\r\nTTTTT\r\n");
            Assert.AreEqual(sml2.Start_Database, -5);
            Assert.AreEqual(sml2.Start_Query, -5);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("C");
            query = new Sequence_ByteArray("G");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 2);
            al = sml2.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml2.Start_Database, 0);
            Assert.AreEqual(sml2.Start_Query, 0);
            Assert.AreEqual(sml2.End_Database, 0);
            Assert.AreEqual(sml2.End_Query, 0);
            Assert.AreEqual(sml2.EndofDynamicExtension, true);
            Assert.AreEqual(sml2.Score, 0.0F);

            database = new Sequence_ByteArray("TTCTTTTTAAAA");
            query = new Sequence_ByteArray("TTTTTTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2.0F, 0.5F), 0);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTCTTTTTAAAA\r\n||.|||||    \r\nTTTTTTTT----\r\n");
            Assert.AreEqual(sml2.Start_Database, -12);
            Assert.AreEqual(sml2.Start_Query, -8);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.AreEqual(sml2.Score, 2.5F);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("TTTTTAAAACCTTT");
            query = new Sequence_ByteArray("TTTTTCCTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 5);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCTTT\r\n|||||\r\nCCTTT\r\n");
            Assert.AreEqual(sml2.Start_Database, -5);
            Assert.AreEqual(sml2.Start_Query, -5);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml2.Score - 5.0F) < 0.001);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("TTTTTAAAACCTTT");
            query = new Sequence_ByteArray("TTTTTCCTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TAAAACCTTT\r\n|    |||||\r\nT----CCTTT\r\n");
            Assert.AreEqual(sml2.Start_Database, -10);
            Assert.AreEqual(sml2.Start_Query, -6);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml2.Score - 4.7F) < 0.0001);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);

            database = new Sequence_String("TTTTTAAAACCTTT");
            query = new Sequence_String("TTTTTCCTTT");
            sml2 = new SmithWatermanGotoh_DynamicBanded_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml2.GetAlignment();
            Assert.AreEqual(al.ToString(), "TAAAACCTTT\r\n|    |||||\r\nT----CCTTT\r\n");
            Assert.AreEqual(sml2.Start_Database, -10);
            Assert.AreEqual(sml2.Start_Query, -6);
            Assert.AreEqual(sml2.End_Database, -1);
            Assert.AreEqual(sml2.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml2.Score - 4.7F) < 0.0001);
            Assert.AreEqual(sml2.EndofDynamicExtension, false);
        }

        [Test]
        public void Test_SmitWatermanGotoh_DynamicBanded_454_3p()
        {
            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            IDynamicBandedDynamicProgramming sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 2), 0);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);
            Assert.AreEqual(sml.ToString(), "\t$\tT\t\n$\t0,00\t-9,00\nT\t-9,00\t2,00\n");


            database = new Sequence_ByteArray("TTT");
            query = new Sequence_ByteArray("TTTATATAT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.End_Database, 3);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("TTTTTTTTTT");
            query = new Sequence_ByteArray("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 10);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_String("TTTTTTTTTT");
            query = new Sequence_String("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.End_Query, 10);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.EndofDynamicExtension, false);





            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("T");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);


            ///
            ///Test for roberts law
            ///
            database = new Sequence_ByteArray("CGATTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTT-CGTCGT\r\n|||||| ||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTCGTCGT\r\n||||||.||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");

            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");

            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 20, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");


            ///
            ///
            ///Unspecific
            ///
            database = new Sequence_ByteArray("ACGTACGTTTTTTTACGTACGTACGAAAAAAATCGTCGTGTA");
            query = new Sequence_ByteArray("ACGAACGTTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTTTTTTTACGTACGTACGAAAAAAA-TCGTCGTGTA\r\n|||.||||||||| |||||||.|||||||||| ||||||.|||\r\nACGAACGTTTTTT-ACGTACGCACGAAAAAAAATCGTCGAGTA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 42);
            Assert.AreEqual(sml.End_Query, 42);
            Assert.IsTrue(Math.Abs(sml.Score - 47.14) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CCCCCC-GGGGGGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| |||||| ||||||| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 41);
            Assert.AreEqual(sml.End_Query, 47);
            Assert.IsTrue(Math.Abs(sml.Score - 106.0F) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CTCCCCCCAGAGAGG--AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||.|||||.|.|.||  ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");


            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_String("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_String("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTTATAAATTAAAA");
            query = new Sequence_ByteArray("CGTCGTTAAATAAAAAAA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTTA--TAAATTAAAA\r\n||||||||  ||||  ||||\r\nCGTCGTTAAATAAA--AAAA\r\n");

            database = new Sequence_ByteArray("CGTCGATTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGAATTTCGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGA-TTTCGTCGT\r\n|||||| |||||||||\r\nCGTCGAATTTCGTCGT\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");

            //
            //Right Cutoff
            //
            database = new Sequence_ByteArray("TTTCCAAAATTTTT");
            query = new Sequence_ByteArray("TTTCCTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 5);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTC\r\n||||\r\nTTTC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 4);
            Assert.AreEqual(sml.End_Query, 4);
            Assert.AreEqual(sml.Score, 4.0F);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("TTTCCAAAATTTTT");
            query = new Sequence_ByteArray("TTTCCTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTCCAAAAT\r\n|||||    |\r\nTTTCC----T\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 10);
            Assert.AreEqual(sml.End_Query, 6);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("CAAAATCGT");
            query = new Sequence_ByteArray("CTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.EndofDynamicExtension, true);

            database = new Sequence_ByteArray("CAAAATCGT");
            query = new Sequence_ByteArray("CTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAAT\r\n|    |\r\nC----T\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 6);
            Assert.AreEqual(sml.End_Query, 2);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("CTAAAATCGT");
            query = new Sequence_ByteArray("CTTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CTAAAAT\r\n||    |\r\nCT----T\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("CTAAAATCGT");
            query = new Sequence_ByteArray("CTTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_3p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "C\r\n|\r\nC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


        }

        [Test]
        public void Test_SmitWatermanGotoh_DynamicBanded_454_5p()
        {
            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            IDynamicBandedDynamicProgramming sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 2), 0);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, -1);
            Assert.AreEqual(sml.Start_Query, -1);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);
            Assert.AreEqual(sml.ToString(), "\t$\tT\t\n$\t0,00\t0,00\nT\t0,00\t2,00\n");


            database = new Sequence_ByteArray("TTT");
            query = new Sequence_ByteArray("TATATATTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(sml.Start_Query, -3);
            Assert.AreEqual(sml.Start_Database, -3);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("TTTTTTTTTT");
            query = new Sequence_ByteArray("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, -10);
            Assert.AreEqual(sml.Start_Database, -10);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_String("TTTTTTTTTT");
            query = new Sequence_String("TTTTTTTTTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2.5F, 5, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TTTTTTTTTT\r\n||||||||||\r\nTTTTTTTTTT\r\n");
            Assert.AreEqual(sml.Start_Query, -10);
            Assert.AreEqual(sml.Start_Database, -10);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("A");
            query = new Sequence_ByteArray("T");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.EndofDynamicExtension, true);


            ///
            ///Test for roberts law
            ///
            database = new Sequence_ByteArray("CGATTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTT-CGTCGT\r\n|||||| ||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTCGTCGT\r\n||||||.||||||\r\nCGATTTACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");

            database = new Sequence_ByteArray("CGATTTTTCGTCGT");
            query = new Sequence_ByteArray("CGATTTACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGATTTTTCGTCGT\r\n|||||| .||||||\r\nCGATTT-ACGTCGT\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(2, 2, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");


            database = new Sequence_ByteArray("CAAAACAACAGAAACAAAACAAAAACACA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(4, 6, 20, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CAAAACAACAG-AAACAAAACAAAAACACA\r\n||||||||||| ||||||||||||||||||\r\nCAAAACAACAGTAAACAAAACAAAAACACA\r\n");

            ///
            ///
            ///Unspecific
            ///
            database = new Sequence_ByteArray("ACGTACGTTTTTTTACGTACGTACGAAAAAAATCGTCGTGTA");
            query = new Sequence_ByteArray("ACGAACGTTTTTTACGTACGCACGAAAAAAAATCGTCGAGTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(9, 1), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTTTTTTTACGTACGTACGAAAAAAA-TCGTCGTGTA\r\n|||.||||||||| |||||||.|||||||||| ||||||.|||\r\nACGAACGTTTTTT-ACGTACGCACGAAAAAAAATCGTCGAGTA\r\n");
            Assert.AreEqual(sml.Start_Database, -42);
            Assert.AreEqual(sml.Start_Query, -42);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml.Score - 47.14) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCCCCCCGGGGGGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CCCCCC-GGGGGGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| |||||| ||||||| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCCCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, -41);
            Assert.AreEqual(sml.Start_Query, -47);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml.Score - 106.0F) < 0.01);

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CTCCCCCCAGAGAGG--AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||.|||||.|.|.||  ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_String("CGTCGTAAAATTTTTCTCCCCAGAGAGGAAAAAAACGTCGT");
            query = new Sequence_String("CGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTAAAA--TTTTT-CT---CCCCAGAGAGG-AAAAAAA-CGTCGT\r\n||||||||||  ||||| ||   ||||.|.|.|| ||||||| ||||||\r\nCGTCGTAAAAAATTTTTTCTACCCCCCGGGGGGGGAAAAAAAACGTCGT\r\n");

            database = new Sequence_ByteArray("CGTCGTTATAAATTAAAA");
            query = new Sequence_ByteArray("CGTCGTTAAATAAAAAAA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTTA--TAAATTAAAA\r\n||||||||  ||||  ||||\r\nCGTCGTTAAATAAA--AAAA\r\n");

            database = new Sequence_ByteArray("CGTCGATTTCGTCGT");
            query = new Sequence_ByteArray("CGTCGAATTTCGTCGT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGA-TTTCGTCGT\r\n|||||| |||||||||\r\nCGTCGAATTTCGTCGT\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");

            database = new Sequence_ByteArray("TATATTAAGTGAAATTTTATATTTAAATTA");
            query = new Sequence_ByteArray("TATATTAAGTGAAATTTTAATATTTAAATTA");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TATATTAAGTGAAATTTTA-TATTTAAATTA\r\n||||||||||||||||||| |||||||||||\r\nTATATTAAGTGAAATTTTAATATTTAAATTA\r\n");


            //
            //Test behaviour
            //
            database = new Sequence_ByteArray("TTTTTAAAACCTTT");
            query = new Sequence_ByteArray("TTTTTCCTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 5);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CCTTT\r\n|||||\r\nCCTTT\r\n");
            Assert.AreEqual(sml.Start_Database, -5);
            Assert.AreEqual(sml.Start_Query, -5);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.IsTrue(Math.Abs(sml.Score - 5.0F) < 0.001);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("TTTTTAAAACCTTT");
            query = new Sequence_ByteArray("TTTTTCCTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 4);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TAAAACCTTT\r\n|    |||||\r\nT----CCTTT\r\n");
            Assert.AreEqual(sml.Start_Database, -10);
            Assert.AreEqual(sml.Start_Query, -6);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);

            database = new Sequence_ByteArray("CGAATTTTAAAA");
            query = new Sequence_ByteArray("CGAATTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 7);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TAAAA\r\n|    \r\nT----\r\n");
            Assert.AreEqual(sml.Start_Database, -5);
            Assert.AreEqual(sml.Start_Query, -1);
            Assert.AreEqual(sml.End_Database, -1);
            Assert.AreEqual(sml.End_Query, -1);
            Assert.AreEqual(sml.EndofDynamicExtension, false);


            database = new Sequence_ByteArray("CGAATTTTAAAA");
            query = new Sequence_ByteArray("CGAATTTT");
            sml = new SmithWatermanGotoh_DynamicBanded_454P_5p(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 1.0F, 0.1F), 8);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.EndofDynamicExtension, true);

        }
    }

    [TestFixture]
    public class Aln_IO
    {
        public static void Start()
        {
            Aln_IO p = new Aln_IO();
            p.Test_PairwiseNucleotideSequenceReader();
            p.Test_PairwiseNucleotideSequenceWriter();

        }
        [Test]
        public void Test_PairwiseNucleotideSequenceReader()
        {
            StreamReader sr = new StreamReader("G:\\Programs\\00TestSequencen\\Nunit_IOAlignment_PairwiseNucleotideSequenceAlignmentReader.txt", Encoding.ASCII);
            PairwiseNucleotideSequenceAlignmentReader nsr = new PairwiseNucleotideSequenceAlignmentReader(sr);
            List<IPairwiseNucleotideSequenceAlignmentContainer> pwa = nsr.GetAllAlignments();
            sr.Close();
            Assert.AreEqual(pwa.Count, 3);
            Assert.AreEqual(pwa[0].Count_SubAlignments, 1);
            Assert.AreEqual(pwa[1].Count_SubAlignments, 4);
            Assert.AreEqual(pwa[2].Count_SubAlignments, 1);
            Assert.AreEqual(pwa[0].DatabaseParent, "data1");
            Assert.AreEqual(pwa[0].QueryParent, "queue1"); //>	10	70	1	42	126,00	100,00	0
            Assert.AreEqual(pwa[0].StartDatabase, 10);
            Assert.AreEqual(pwa[0].EndDatabase, 70);
            Assert.AreEqual(pwa[0].StartQuery, 1);
            Assert.AreEqual(pwa[0].EndQuery, 42);
            Assert.AreEqual(pwa[0].Score, 126.0F);
            Assert.AreEqual(pwa[0].LengthDatabaseParent, 2500);
            Assert.AreEqual(pwa[0].LengthQueryParent, 250);
            Assert.AreEqual(pwa[1].DatabaseParent, "data2");
            Assert.AreEqual(pwa[1].QueryParent, "queue2");//10	90	1	108	318,00	100,00	1
            Assert.AreEqual(pwa[1].StartDatabase, 10);
            Assert.AreEqual(pwa[1].StartQuery, 1);//500	600	1	73	219,00	100,00	0
            Assert.AreEqual(pwa[1].EndDatabase, 600);
            Assert.AreEqual(pwa[1].EndQuery, 73);
            Assert.AreEqual(pwa[1][1].StartDatabase, 120);//120	200	1	73	219,00	100,00	0
            Assert.AreEqual(pwa[1][1].EndDatabase, 200);
            Assert.AreEqual(pwa[1][1].Score, 219.0F);
            Assert.AreEqual(pwa[1][1].PlusPlusStrand, true);
            Assert.AreEqual(pwa[1][0].Count_LongGaps_Database(1), 1);
            Assert.AreEqual(pwa[1].Count_LongGaps_Query(10), 3);
            Assert.AreEqual(pwa[1][3].Alignment.DatabaseSequence.ToString(), "CGAGAGGTGAAATTCTTGGACCGTCGTAAGACTAACTTAAGCGAAAGCATTTGCCAAAGATGTTTTCATTGAT");
            Assert.AreEqual(pwa[1][3].Alignment.SimilaritySequence.ToString(), "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||");
            Assert.AreEqual(pwa[1][3].Alignment.QuerySequence.ToString(), "CGAGAGGTGAAATTCTTGGACCGTCGTAAGACTAACTTAAGCGAAAGCATTTGCCAAAGATGTTTTCATTAAT");
            Assert.AreEqual(pwa[1].LengthDatabaseParent, 3500);
            Assert.AreEqual(pwa[1].LengthQueryParent, 350);

            Assert.AreEqual(pwa[2].Alignment.DatabaseSequence.ToString(), "GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT");
            Assert.AreEqual(pwa[2].Alignment.QuerySequence.ToString(), "GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT");
            Assert.AreEqual(pwa[2].LengthDatabaseParent, 4500);
            Assert.AreEqual(pwa[2].LengthQueryParent, 450);


            Assert.AreEqual(pwa[0].PlusPlusStrand, true);
            Assert.AreEqual(pwa[1].PlusPlusStrand, true);
            Assert.AreEqual(pwa[2].PlusPlusStrand, false);
            Assert.AreEqual(pwa[2].DatabaseParent, "data3");
            Assert.AreEqual(pwa[2].QueryParent, "queue3");

            //
            //ROUND2
            //
            pwa.Clear();
            sr = new StreamReader("G:\\Programs\\00TestSequencen\\Nunit_IOAlignment_PairwiseNucleotideSequenceAlignmentReader.txt", Encoding.ASCII);
            nsr = new PairwiseNucleotideSequenceAlignmentReader(sr, 1);
            int count = 0;
            IPairwiseNucleotideSequenceAlignmentContainer temp;
            while ((temp = nsr.GetNextAlignment())!=null)
            {
                pwa.Add(temp);
                count++;
            }
            Assert.AreEqual(pwa.Count, 3);
            Assert.AreEqual(count, 3);
            Assert.AreEqual(pwa[0].DatabaseParent, "data1");
            Assert.AreEqual(pwa[1].DatabaseParent, "data2");
            Assert.AreEqual(pwa[2].DatabaseParent, "data3");




        }

        [Test]
        public void Test_PairwiseNucleotideSequenceWriter()
        {

            PairwiseNucleotideSequenceAlignment pw1 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAATT"), new Sequence_ByteArray("AATTT")), "data1", "queue1", 23, 123, 93, 192);
            pw1.PlusPlusStrand = false;
            pw1.Score = 12;

            PairwiseNucleotideSequenceAlignment pw2 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCC"), new Sequence_ByteArray("GGG")), "data2", "queue2", 12, 123, 112, 192);
            pw2.PlusPlusStrand = true;
            pw2.Score = 15;

            PairwiseNucleotideSequenceAlignment pw3 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT"), new Sequence_ByteArray("GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT")), "data2", "queue2", 212, 56, 252, 156);
            pw3.PlusPlusStrand = true;
            pw3.Score = 40;

            PairwiseNucleotideSequenceAlignment pw4 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAT"), new Sequence_ByteArray("AAC")), "data2", "queue2", 312, 123, 353, 345);
            pw4.PlusPlusStrand = true;
            pw4.Score = 11;

            PairwiseNucleotideSequenceAlignment pw5 = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ACGACG"), new Sequence_ByteArray("ACGACG")), "data3", "queue3", 312, 123, 353, 345);
            pw5.PlusPlusStrand = true;
            pw5.Score = 11;
            pw5.LengthDatabaseParent = 2500;
            pw5.LengthQueryParent = 465;
            List<IPairwiseNucleotideSequenceAlignmentContainer> pwc = new List<IPairwiseNucleotideSequenceAlignmentContainer>();
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            temp.Add(pw2);
            temp.Add(pw3);
            temp.Add(pw4);
            pwc.Add(pw1);
            pwc.Add(new CompositePairwiseNucleotideSequenceAlignment(temp));
            pwc.Add(pw5);

            StreamWriter sw = new StreamWriter("G:\\Programs\\00TestSequencen\\Nunit_IOAlignment_PairwiseNucleotideSequenceAlignmentWriter.txt", false, Encoding.ASCII);
            PairwiseNucleotideSequenceAlignmentWriter nswr = new PairwiseNucleotideSequenceAlignmentWriter(sw);
            nswr.WriteAlignments(pwc);
            sw.Close();

            StreamReader sr = new StreamReader("G:\\Programs\\00TestSequencen\\Nunit_IOAlignment_PairwiseNucleotideSequenceAlignmentWriter.txt", Encoding.ASCII);
            PairwiseNucleotideSequenceAlignmentReader nsr = new PairwiseNucleotideSequenceAlignmentReader(sr);
            List<IPairwiseNucleotideSequenceAlignmentContainer> pwa = nsr.GetAllAlignments();
            sr.Close();

            Assert.AreEqual(pwa.Count, 3);
            Assert.AreEqual(pwa[0].Count_SubAlignments, 1);
            Assert.AreEqual(pwa[1].Count_SubAlignments, 3);
            Assert.AreEqual(pwa[2].Count_SubAlignments, 1);
            Assert.AreEqual(pwa[0].DatabaseParent, "data1");
            Assert.AreEqual(pwa[1].DatabaseParent, "data2");
            Assert.AreEqual(pwa[2].DatabaseParent, "data3");
            Assert.AreEqual(pwa[0].QueryParent, "queue1");
            Assert.AreEqual(pwa[1].QueryParent, "queue2");
            Assert.AreEqual(pwa[2].QueryParent, "queue3");
            Assert.AreEqual(pwa[0].PlusPlusStrand, false);
            Assert.AreEqual(pwa[1].PlusPlusStrand, true);
            Assert.AreEqual(pwa[2].PlusPlusStrand, true);
            Assert.AreEqual(pwa[0].Score, 12);
            Assert.AreEqual(pwa[1].Score, 66);
            Assert.AreEqual(pwa[2].Score, 11);

            Assert.AreEqual(pwa[0].Alignment.DatabaseSequence.ToString(), "AAATT");
            Assert.AreEqual(pwa[0].Alignment.QuerySequence.ToString(), "AATTT");
            Assert.AreEqual(pwa[1][1].Alignment.DatabaseSequence.ToString(), "GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT");
            Assert.AreEqual(pwa[1][1].Alignment.QuerySequence.ToString(), "GCCGGCCAACGTGGCGGTCAACGGGTCGTCCAACGCGGAATTCCACGCGCTGTTGCAGGCGACGAGGATGACGAGGACGACGTGAACATCGACGATATCTAAAAGGAATCAAAAAGAGAACAGCTACCCTAATCCAACCAGTAATATTTTAAAAGGGGGAAAATCAAGCGAGCAAGCATAATCTACAATTCAAACTGAAATATAAATTAGCGAGTTAAATCAACAT");
            Assert.AreEqual(pwa[1][1].StartDatabase, 212);//	212	252	56	156	40,00	100,00	0
            Assert.AreEqual(pwa[1][1].StartQuery, 56);
            Assert.AreEqual(pwa[1][1].EndDatabase, 252);
            Assert.AreEqual(pwa[1][1].EndQuery, 156);
            Assert.AreEqual(pwa[1][1].PlusPlusStrand, true);
            Assert.AreEqual(pwa[1].LengthDatabaseParent, null);
            Assert.AreEqual(pwa[1].LengthQueryParent, null);


            Assert.AreEqual(pwa[1][0].Alignment.DatabaseSequence.ToString(), "CCC");
            Assert.AreEqual(pwa[1][0].Alignment.QuerySequence.ToString(), "GGG");
            Assert.AreEqual(pwa[1][2].Alignment.DatabaseSequence.ToString(), "AAT");
            Assert.AreEqual(pwa[1][2].Alignment.QuerySequence.ToString(), "AAC");
            Assert.AreEqual(pwa[1][2].Score, 11);//312	353	123	345	11,00	66,67	0
            Assert.AreEqual(pwa[1][2].StartDatabase, 312);
            Assert.AreEqual(pwa[1][2].StartQuery, 123);
            Assert.AreEqual(pwa[1][2].EndDatabase, 353);
            Assert.AreEqual(pwa[1][2].EndQuery, 345);

            Assert.AreEqual(pwa[2].Alignment.DatabaseSequence.ToString(), "ACGACG");
            Assert.AreEqual(pwa[2].Alignment.QuerySequence.ToString(), "ACGACG");
            Assert.AreEqual(pwa[2].StartDatabase, 312);
            Assert.AreEqual(pwa[2].StartQuery, 123);
            Assert.AreEqual(pwa[2].EndDatabase, 353); //	312	353	123	345	11,00	100,00	0
            Assert.AreEqual(pwa[2].EndQuery, 345);

            Assert.AreEqual(pwa[2].LengthDatabaseParent, 2500);
            Assert.AreEqual(pwa[2].LengthQueryParent, 465);



        }
    }

    [TestFixture]
    public class Aln_Blast
    {
        public static void Start()
        {
            Aln_Blast p = new Aln_Blast();
            p.Test_BlastHashProcessor();
            p.Test_BlastPostProcessor_PartialAlginmentAggregator();
            p.Test_BlastSeedProcessor_BestCumulativeDiagonal();
            p.Test_BlastSeedProcessor_BestDiagonal();
            p.Test_IntronPolisher_KeepBestOverlap();
            p.Test_Anchored454SmithWatermanGotoh();
            p.Test_AnchoredDynamicBandedSmithWatermanGotoh();
            p.Test_AnchoredDynamicBanded454SmithWatermanGotoh();
        }


        [Test]
        public void Test_BlastHashProcessor()
        {
            //
            //Non overlapping blast hash processor
            //
            BlastHashProcessor_NonOverlapping bh = new BlastHashProcessor_NonOverlapping(3, 100);
            NucleotideSequence ns1 = new NucleotideSequence("fritzi", new Sequence_ByteArray("TTATTACCGCCGACGAGCA")); //length=19
            NucleotideSequence ns2 = new NucleotideSequence("franz", new Sequence_ByteArray("CGCGTCGTCGTCGTCACCCCC"));//length=21


            DNADictionary<List<BlastHashSeed>> dict = bh.GetBlastHashtable(ns1);
            List<BlastHashSeed> bhs = dict[new Sequence_ByteArray("TTA")];

            Assert.AreEqual(bhs.Count, 2);
            Assert.AreEqual(bhs[0].DatabaseID, 0);
            Assert.AreEqual(bhs[0].DatabasePosition, 0);
            Assert.AreEqual(bhs[1].DatabaseID, 0);
            Assert.AreEqual(bhs[1].DatabasePosition, 3);

            //Test removal of low complexity
            ns1 = new NucleotideSequence("fritzi", new Sequence_ByteArray("TTATTATTATTATTACCGCCGACGAGCA")); //length=19
            bh = new BlastHashProcessor_NonOverlapping(3,2);
            dict = bh.GetBlastHashtable(ns1);
            bhs = dict[new Sequence_ByteArray("TTA")];
            Assert.AreEqual(bhs, null);
            //Removal of low complexity regions work


            ///
            ///Test the list version
            ///
            List<NucleotideSequence> list = new List<NucleotideSequence>();
            list.Add(ns1);
            list.Add(ns2);
            bh = new BlastHashProcessor_NonOverlapping(3, 100);
            dict = bh.GetBlastHashtable(list);
            bhs = dict[new Sequence_ByteArray("TTA")];
            Assert.AreEqual(bhs.Count, 5);
            Assert.AreEqual(bhs[0].DatabaseID, 0);
            Assert.AreEqual(bhs[0].DatabasePosition, 0);
            Assert.AreEqual(bhs[1].DatabaseID, 0);
            Assert.AreEqual(bhs[1].DatabasePosition, 3);

            bhs = dict[new Sequence_ByteArray("GTC")];
            Assert.AreEqual(bhs.Count, 4);
            Assert.AreEqual(bhs[0].DatabaseID, 1);
            Assert.AreEqual(bhs[0].DatabasePosition, 3);
            Assert.AreEqual(bhs[1].DatabaseID, 1);
            Assert.AreEqual(bhs[1].DatabasePosition, 6);
            Assert.AreEqual(bhs[2].DatabaseID, 1);
            Assert.AreEqual(bhs[2].DatabasePosition, 9);
            Assert.AreEqual(bhs[3].DatabaseID, 1);
            Assert.AreEqual(bhs[3].DatabasePosition, 12);





        }

        [Test]
        public void Test_BlastSeedProcessor_BestDiagonal()
        {

            BlastSeedProcessor_BestDiagonal sP;
            List<BlastSeed> seeds;
            List<BlastSeed> results;

            sP = new BlastSeedProcessor_BestDiagonal(2, 2);
            seeds = new List<BlastSeed>();
            seeds.Add(new BlastSeed(12, 3));
            seeds.Add(new BlastSeed(15, 5));
            results = sP.GetSeeds(seeds);

            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseID, 0);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].QueryStartPosition, 3);
            Assert.AreEqual(results[0].Offset, 9.0F);
            Assert.AreEqual(results[0].Score, 2);


            sP = new BlastSeedProcessor_BestDiagonal(2, 2, 2);
            seeds = new List<BlastSeed>();
            seeds.Add(new BlastSeed(2, 15, 10, 5));
            seeds.Add(new BlastSeed(2, 12, 10, 3));
            results = sP.GetSeeds(seeds);

            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseID, 2);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].QueryStartPosition, 3);
            Assert.AreEqual(results[0].Offset, 9.0F);
            Assert.AreEqual(results[0].Score, 2);


            seeds.Add((new BlastSeed(2, 17, 10, 6)));
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].Score, 3);


            //Add seed with a offset which is to large to be considered
            seeds.Add((new BlastSeed(2, 20, 4)));
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].Score, 3);

            seeds.Add(new BlastSeed(4, 15, 5));
            seeds.Add(new BlastSeed(4, 12, 3));
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseID, 2);
            Assert.AreEqual(results[0].Score, 3);
            Assert.AreEqual(results[1].Score, 2);
            Assert.AreEqual(results[1].DatabaseID, 4);

            ///TEst with equal score, will it give back more than one sequence
            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 2, 0);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseID, 2);

            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 2, 1);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);

            seeds.Add(new BlastSeed(4, 16, 6));
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseID, 4);
            Assert.AreEqual(results[0].Score, 3);
            Assert.AreEqual(results[1].Score, 3);
            Assert.AreEqual(results[1].DatabaseID, 2);


            seeds.Add(new BlastSeed(5, 15, 10, 5));
            seeds.Add(new BlastSeed(5, 12, 10, 3));
            seeds.Add(new BlastSeed(5, 16, 11, 6));
            seeds.Add(new BlastSeed(6, 12, 10, 1));

            ///Check if it returns all possible seeds even when it will return 3 instead of 1
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 3);
            Assert.AreEqual(results[0].DatabaseID, 5);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].QueryStartPosition, 3);

            ///Test if it returns all seeds fullfilling the minimum requirements
            sP = new BlastSeedProcessor_BestDiagonal(2, true, 5);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 3);


            //Test if even works for single seeds
            sP = new BlastSeedProcessor_BestDiagonal(1, true, 5);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 4);


            //New Round
            seeds.Clear();
            seeds.Add(new BlastSeed(5, 15, 10, 5));
            seeds.Add(new BlastSeed(5, 16, 10, 6));
            seeds.Add(new BlastSeed(5, 17, 10, 7));
            seeds.Add(new BlastSeed(5, 18, 10, 8));
            seeds.Add(new BlastSeed(6, 16, 10, 6));
            seeds.Add(new BlastSeed(6, 12, 10, 2));
            seeds.Add(new BlastSeed(7, 12, 10, 2));
            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 5, 1);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 5, 0);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 5, 2);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            seeds.Add(new BlastSeed(7, 13, 10, 3));
            sP = new BlastSeedProcessor_BestDiagonal(2, 1, 5, 2);
            results = sP.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 3);












        }

        [Test]
        public void Test_BlastSeedProcessor_BestCumulativeDiagonal()
        {
            //Initialize working variables
            List<BlastSeed> seeds = new List<BlastSeed>();
            List<BlastSeed> results = new List<BlastSeed>();
            IBlastHashtableConstructor hashConst = new BlastHashProcessor_NonOverlapping(10, 100);
            BlastSeedProcessor_BestCumulativeDiagonal seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(3, 10, hashConst, null);



            seeds.Add(new BlastSeed(2, 20, 120));
            seeds.Add(new BlastSeed(2, 12, 112));
            seeds.Add(new BlastSeed(2, 15, 115));

            seeds.Add(new BlastSeed(3, 27, 126));
            seeds.Add(new BlastSeed(3, 24, 123));

            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].Score, 3);
            Assert.AreEqual(results[0].DatabaseStartPosition, 12);
            Assert.AreEqual(results[0].DatabaseEndPosition, 41);
            Assert.AreEqual(results[0].QueryStartPosition, 112);
            Assert.AreEqual(results[0].QueryEndPosition, 141);
            Assert.AreEqual(results[0].Offset, -100);

            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 10, hashConst, null);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);


            seeds.Add(new BlastSeed(3, 50, 223));
            seeds.Add(new BlastSeed(3, 55, 228));
            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, null, 5, 0);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseID, 3);
            Assert.AreEqual(results[1].DatabaseID, 3);
            Assert.AreEqual(results[0].DatabaseStartPosition, 24);
            Assert.AreEqual(results[0].DatabaseEndPosition, 43);
            Assert.AreEqual(results[1].DatabaseStartPosition, 50);
            Assert.AreEqual(results[1].DatabaseEndPosition, 69);

            seeds.Clear();
            seeds.Add(new BlastSeed(2, 40, 220));
            seeds.Add(new BlastSeed(2, 20, 120));
            seeds.Add(new BlastSeed(2, 50, 230));
            seeds.Add(new BlastSeed(2, 12, 112));
            seeds.Add(new BlastSeed(2, 15, 115));
            seeds.Add(new BlastSeed(2, 60, 240));

            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(3, 1, hashConst, null);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[1].DatabaseStartPosition, 12);
            Assert.AreEqual(results[1].DatabaseEndPosition, 41);
            Assert.AreEqual(results[0].DatabaseStartPosition, 40);
            Assert.AreEqual(results[0].DatabaseEndPosition, 69);

            ///Test the maximum distance
            seeds.Clear();
            seeds.Add(new BlastSeed(3, 1, 100));
            seeds.Add(new BlastSeed(3, 11, 111));
            seeds.Add(new BlastSeed(3, 2020, 121));
            seeds.Add(new BlastSeed(3, 2030, 131));
            seeds.Add(new BlastSeed(3, 2040, 141));

            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2001);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            Assert.AreEqual(results[0].DatabaseStartPosition, 1);
            Assert.AreEqual(results[1].DatabaseStartPosition, 2020);


            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2000, 5, 0);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            Assert.AreEqual(results[0].DatabaseStartPosition, 2020);
            //Explanation: in the first scenary two diagonals are found which will be aggregated into a single cumulative diagonal.
            //This cumulative diagonal is split into two sepparate alignments in the end to allow two independent local alignments.
            //In the second scenary the two diagonals are not aggregated and only the single best will be returned which is the diagonal encompassing three seeds.


            //
            //NEW Round
            //
            seeds.Clear();
            seeds.Add(new BlastSeed(5, 40, 240));
            seeds.Add(new BlastSeed(5, 50, 250));
            seeds.Add(new BlastSeed(6, 250, 50));
            seeds.Add(new BlastSeed(2, 20, 120));
            seeds.Add(new BlastSeed(2, 12, 112));
            seeds.Add(new BlastSeed(2, 15, 115));
            seeds.Add(new BlastSeed(2, 40, 140));

            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2000, 5, 0);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2000, 5, 1);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 1);
            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2000, 5, 2);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 2);
            seeds.Add(new BlastSeed(6, 260, 60));
            seedProcessor = new BlastSeedProcessor_BestCumulativeDiagonal(2, 1, hashConst, 2000, 5, 2);
            results = seedProcessor.GetSeeds(seeds);
            Assert.AreEqual(results.Count, 3);





        }

        [Test]
        public void Test_IntronPolisher_KeepBestOverlap()
        {
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            List<PairwiseNucleotideSequenceAlignment> res = new List<PairwiseNucleotideSequenceAlignment>();
            IntronBoundaryPolishing_KeepBestOverlap intronPolisher;
            PairwiseNucleotideSequenceAlignment pw;

            intronPolisher = new IntronBoundaryPolishing_KeepBestOverlap(SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(5.0F, 3.0F, 11.0F, 1.5F), 1);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCA"), new Sequence_ByteArray("CCCA")), "data1", "query1", 1, 3, 4, 6);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 4, 6, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCA\r\n||||\r\nCCCA\r\n");
            Assert.AreEqual(res[0].PlusPlusStrand, true);
            Assert.AreEqual(res[0].DatabaseParent, "data1");
            Assert.AreEqual(res[0].QueryParent, "query1");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].EndDatabase, 4);
            Assert.AreEqual(res[0].EndQuery, 6);

            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(res[1].PlusPlusStrand, true);
            Assert.AreEqual(res[1].DatabaseParent, "data1");
            Assert.AreEqual(res[1].QueryParent, "query1");
            Assert.AreEqual(res[1].StartDatabase, 5);
            Assert.AreEqual(res[1].StartQuery, 7);
            Assert.AreEqual(res[1].EndDatabase, 7);
            Assert.AreEqual(res[1].EndQuery, 9);


            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCAAA"), new Sequence_ByteArray("CCCATA")), "data1", "query1", 1, 3, 6, 8);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAATTT"), new Sequence_ByteArray("AAATTT")), "data1", "query1", 4, 6, 9, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCC\r\n|||\r\nCCC\r\n");
            Assert.AreEqual(res[0].PlusPlusStrand, true);
            Assert.AreEqual(res[0].DatabaseParent, "data1");
            Assert.AreEqual(res[0].QueryParent, "query1");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].EndDatabase, 3);
            Assert.AreEqual(res[0].EndQuery, 5);

            Assert.AreEqual(res[1].Alignment.ToString(), "AAATTT\r\n||||||\r\nAAATTT\r\n");
            Assert.AreEqual(res[1].PlusPlusStrand, true);
            Assert.AreEqual(res[1].DatabaseParent, "data1");
            Assert.AreEqual(res[1].QueryParent, "query1");
            Assert.AreEqual(res[1].StartDatabase, 4);
            Assert.AreEqual(res[1].StartQuery, 6);
            Assert.AreEqual(res[1].EndDatabase, 9);
            Assert.AreEqual(res[1].EndQuery, 11);


            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAA"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AATATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 4, 6, 9, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCATAA\r\n||||.||\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AACGATCGATAATTCC"), new Sequence_ByteArray("AACGATCGATAAATCC")), "data1", "query1", 1, 3, 16, 18);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AATTCCTTCGTACGTACGTA"), new Sequence_ByteArray("AATTCCTTCGTACGTACGTA")), "data1", "query1", 11, 13, 30, 32);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);
            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "AACGATCGAT\r\n||||||||||\r\nAACGATCGAT\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "AATTCCTTCGTACGTACGTA\r\n||||||||||||||||||||\r\nAATTCCTTCGTACGTACGTA\r\n");


            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAA"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AATAATTT"), new Sequence_ByteArray("A-AAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCATAA\r\n||||.||\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 19);
            Assert.AreEqual(res[0].EndDatabase, 7);
            Assert.AreEqual(res[1].EndDatabase, 21);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAAA"), new Sequence_ByteArray("CCCA-AAA")), "data1", "query1", 1, 3, 8, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 14, 6, 20, 12);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCC\r\n|||\r\nCCC\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "ATAATTT\r\n|.|||||\r\nAAAATTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 14);
            Assert.AreEqual(res[0].EndDatabase, 3);
            Assert.AreEqual(res[1].EndDatabase, 20);



            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAA"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AA-ATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCATAA\r\n||||.||\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 17);
            Assert.AreEqual(res[0].EndDatabase, 7);
            Assert.AreEqual(res[1].EndDatabase, 21);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAA"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AAAAAATTT"), new Sequence_ByteArray("AA--AATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCATAA\r\n||||.||\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 20);
            Assert.AreEqual(res[0].EndDatabase, 7);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 10);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATAA"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("A-AAAATTT"), new Sequence_ByteArray("AA--AATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCATAA\r\n||||.||\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "TTT\r\n|||\r\nTTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 19);
            Assert.AreEqual(res[0].EndDatabase, 7);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 10);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCAA-A"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 6, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCC\r\n|||\r\nCCC\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "ATAATTT\r\n|.|||||\r\nAAAATTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 14);
            Assert.AreEqual(res[0].EndDatabase, 3);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 6);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CC-AA-A"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 5, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CC\r\n||\r\nCC\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "ATAATTT\r\n|.|||||\r\nAAAATTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 14);
            Assert.AreEqual(res[0].EndDatabase, 2);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 6);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCAA-A"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 6, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("AAAA--T")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCAA-A\r\n||||| |\r\nCCCAAAA\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 20);
            Assert.AreEqual(res[0].EndDatabase, 6);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 10);

            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCGAA-A"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 6, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("AAAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CC\r\n||\r\nCC\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "ATAATTT\r\n|.|||||\r\nAAAATTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 14);
            Assert.AreEqual(res[0].EndDatabase, 2);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 6);

            //
            //Round 3 - one partial alignment will be entirely removed
            //
            temp.Clear();
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCGAA-A"), new Sequence_ByteArray("CCCAAAA")), "data1", "query1", 1, 3, 6, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("AT"), new Sequence_ByteArray("AT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATAATTT"), new Sequence_ByteArray("ATAATTT")), "data1", "query1", 14, 6, 21, 11);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            pw.LengthDatabaseParent = 230;
            pw.LengthQueryParent = 110;
            temp.Add(pw);
            res = intronPolisher.GetPolishedAlignments(temp);

            Assert.AreEqual(res.Count, 2);
            Assert.AreEqual(res[0].Alignment.ToString(), "CC\r\n||\r\nCC\r\n");
            Assert.AreEqual(res[1].Alignment.ToString(), "ATAATTT\r\n|||||||\r\nATAATTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[1].StartDatabase, 14);
            Assert.AreEqual(res[0].EndDatabase, 2);
            Assert.AreEqual(res[1].EndDatabase, 21);
            Assert.AreEqual(res[1].StartQuery, 6);

        }

        [Test]
        public void Test_BlastPostProcessor_PartialAlginmentAggregator()
        {
            List<PairwiseNucleotideSequenceAlignment> temp = new List<PairwiseNucleotideSequenceAlignment>();
            List<IPairwiseAlignmentContainer> res = new List<IPairwiseAlignmentContainer>();
            SubstitutionMatrix substitutionMatrix = SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F);
            BlastPostProcessor_PartialAlignmentAggregator postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 0);
            PairwiseNucleotideSequenceAlignment pw;


            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCA"), new Sequence_ByteArray("CCCA")), "data1", "query1", 1, 3, 4, 6);
            pw.PlusPlusStrand = true;
            pw.Score = 11;
            temp.Add(pw);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data2", "query1", 4, 6, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].DatabaseParent, "data1");
            Assert.AreEqual(res[0].Score, 11);


            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 2);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data2", "query1", 4, 6, 7, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 12;
            temp.Add(pw);
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 2);
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 2);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 3);

            //
            //ROUND 2
            //
            temp.Clear();
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCA"), new Sequence_ByteArray("CCCA")), "data1", "query1", 1, 3, 4, 6);
            pw.PlusPlusStrand = true;

            pw.Score = 16;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 50, 6, 57, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 58, 10, 98, 34);
            pw.PlusPlusStrand = false;
            pw.Score = 20;
            temp.Add(pw);

            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 0);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCANNNNNNNNNNTTT\r\n||||          |||\r\nCCCA----------TTT\r\n");
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].EndDatabase, 57);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].EndQuery, 9);
            Assert.AreEqual(res[0].PlusPlusStrand, true);
            Assert.AreEqual(res[0].Count_SubAlignments, 2);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 105, 40, 110, 49);
            pw.PlusPlusStrand = true;
            pw.Score = 20;
            temp.Add(pw);
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 0);
            res = postProcessor.GetAlignments(temp);

            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].EndDatabase, 57);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 105, 35, 111, 39);
            pw.PlusPlusStrand = false;
            pw.Score = 30;
            temp.Add(pw);

            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 0);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 58);
            Assert.AreEqual(res[0].EndDatabase, 111);


            //
            //Round3
            //
            temp.Clear();
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCATA"), new Sequence_ByteArray("CCCAAA")), "data1", "query1", 1, 3, 6, 8);
            pw.PlusPlusStrand = true;
            pw.Score = 16;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATATTT"), new Sequence_ByteArray("ATATTT")), "data1", "query1", 50, 6, 57, 13);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 0);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].Count_SubAlignments, 2);
            Assert.AreEqual(res[0].EndDatabase, 57);
            Assert.AreEqual(res[0].EndQuery, 13);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCNNNNNNNNNNATATTT\r\n|||          ||||||\r\nCCC----------ATATTT\r\n");

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("TTCCC"), new Sequence_ByteArray("TTCCC")), "data1", "query1", 100, 12, 107, 16);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);
            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].Count_SubAlignments, 3);
            Assert.AreEqual(res[0].EndDatabase, 107);
            Assert.AreEqual(res[0].EndQuery, 16);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCNNNNNNNNNNATATTTNNNNNNNNNNCCC\r\n|||          ||||||          |||\r\nCCC----------ATATTT----------CCC\r\n");

            //
            //Round4
            //
            temp.Clear();
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCA"), new Sequence_ByteArray("CCCA")), "data1", "query1", 1, 3, 4, 6);
            pw.PlusPlusStrand = true;
            pw.Score = 16;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATATTT"), new Sequence_ByteArray("ATATTT")), "data1", "query1", 40, 50, 57, 13);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATATTT"), new Sequence_ByteArray("ATATTT")), "data1", "query1", 50, 60, 57, 13);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 101, 6, 104, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].Count_SubAlignments, 2);
            Assert.AreEqual(res[0].EndDatabase, 104);
            Assert.AreEqual(res[0].EndQuery, 9);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCANNNNNNNNNNTTT\r\n||||          |||\r\nCCCA----------TTT\r\n");

            //
            //Round5
            //
            temp.Clear();
            postProcessor = new BlastPostProcessor_PartialAlignmentAggregator(1, substitutionMatrix, null, 1);
            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("CCCA"), new Sequence_ByteArray("CCCA")), "data1", "query1", 1, 3, 4, 6);
            pw.PlusPlusStrand = true;
            pw.Score = 16;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATATTT"), new Sequence_ByteArray("ATATTT")), "data1", "query1", 40, 50, 57, 13);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATATTT"), new Sequence_ByteArray("ATATTT")), "data1", "query1", 50, 10, 57, 13);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            pw = new PairwiseNucleotideSequenceAlignment(new PairwiseAlignment(new Sequence_ByteArray("ATTT"), new Sequence_ByteArray("ATTT")), "data1", "query1", 101, 6, 104, 9);
            pw.PlusPlusStrand = true;
            pw.Score = 10;
            temp.Add(pw);

            res = postProcessor.GetAlignments(temp);
            Assert.AreEqual(res.Count, 1);
            Assert.AreEqual(res[0].StartDatabase, 1);
            Assert.AreEqual(res[0].StartQuery, 3);
            Assert.AreEqual(res[0].Count_SubAlignments, 2);
            Assert.AreEqual(res[0].EndDatabase, 57);
            Assert.AreEqual(res[0].EndQuery, 13);
            Assert.AreEqual(res[0].Alignment.ToString(), "CCCANNNNNNNNNNATATTT\r\n||||          ||||||\r\nCCCANNN-------ATATTT\r\n");
        }

        [Test]
        public void Test_Anchored454SmithWatermanGotoh()
        {
            Anchored454SmithWatermanGotoh prototype = new Anchored454SmithWatermanGotoh(0);
            IAnchoredDynamicProgramming sml;

            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            sml = prototype.NewDynamicProgramming(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 1, 1, true);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.Score, 3.0F);


            database = new Sequence_ByteArray("AAAAAGAAAAA");
            query = new Sequence_ByteArray("AGA");
            sml = prototype.NewDynamicProgramming(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 6, 2, true);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AGA\r\n|||\r\nAGA\r\n");
            Assert.AreEqual(sml.Start_Database, 5);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 7);
            Assert.AreEqual(sml.End_Query, 3);
            Assert.AreEqual(sml.Score, 9.0F);



            prototype = new Anchored454SmithWatermanGotoh(75);
            database = new Sequence_ByteArray("AAAAACGTAAAAA");
            query = new Sequence_ByteArray("CGTA");
            sml = prototype.NewDynamicProgramming(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 6, 4, true);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTA\r\n||||\r\nCGTA\r\n");
            Assert.AreEqual(sml.Start_Database, 6);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 9);
            Assert.AreEqual(sml.End_Query, 4);
            Assert.AreEqual(sml.Score, 12.0F);

            prototype = new Anchored454SmithWatermanGotoh(0);
            database = new Sequence_ByteArray("CGTCGTGAAACAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTGAAAACAAAACGTCGT");
            sml = prototype.NewDynamicProgramming(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 1, 1, true);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTGAAA-CAAA-CGTCGT\r\n|||||||||| |||| ||||||\r\nCGTCGTGAAAACAAAACGTCGT\r\n");

            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 20);
            Assert.AreEqual(sml.End_Query, 22);


            database = new Sequence_ByteArray("CGTCGTGAAACAAACGTCGT");
            query = new Sequence_ByteArray("CGTCGTGAAAACAAAACGTCGT");
            sml = prototype.NewDynamicProgramming(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 1.5F), 1, 1, false);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "CGTCGTG-AAAC-AAACGTCGT\r\n||||||| |||| |||||||||\r\nCGTCGTGAAAACAAAACGTCGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 20);
            Assert.AreEqual(sml.End_Query, 22);

        }

        [Test]
        public void Test_AnchoredDynamicBandedSmithWatermanGotoh()
        {
            ISequenceContainer database = new Sequence_ByteArray("A");
            ISequenceContainer query = new Sequence_ByteArray("A");
            AnchoredDynamicBandedSmithWatermanGotoh sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 5, 1), 1, 1, 3, 1, 0);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A\r\n|\r\nA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 1);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            database = new Sequence_ByteArray("TTTTTTTTTTAAAACCGGGGGGGTTTTTTTTTTT"); //Start 14
            query = new Sequence_ByteArray("AAAAGGGGGGG"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1), 14, 4, 10, 2, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAAACCGGGGGGG\r\n||||  |||||||\r\nAAAA--GGGGGGG\r\n");
            Assert.AreEqual(sml.Start_Database, 11);
            Assert.AreEqual(sml.Score, 8);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 23);
            Assert.AreEqual(sml.End_Query, 11);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            database = new Sequence_ByteArray("ACGTTTACCGGGTT"); //Start 14
            query = new Sequence_ByteArray("ACCGTTTACCCGGTT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1), 1, 1, 6, 1, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "A-CGTTTACCGGGTT\r\n| ||||||||.||||\r\nACCGTTTACCCGGTT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 10);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 14);
            Assert.AreEqual(sml.End_Query, 15);
            Assert.AreEqual(sml.UsedDynamicBanding, true);



            database = new Sequence_ByteArray("TCAGTGCCAAACCTACTCCTCAGGCAAATATTGAGCTCTTCTAGGAACTCTTTAACAGCAGTAGTGCAATCACCAGATATACACTCTGCAATGTAGAACGACACTTTTTCTAAGGTAGTTGTGAATCGTAAAGGTCACCTTTTCACTATCTTGAATTCACCAAAAGATCAAAGTAACCAAAGTATCACCAATGACGTGCC"); //Start 14
            query = new Sequence_ByteArray("TCAGTGACAAACCTACGGCTCAGGCAAATATAGAGCTTTACTAGGAACTCTATAACAGCAGCATTGCAATCACCAGATATACACTCTGCTATGTAGAACGACACTTTTTTTAAAGTAGCTGTGAATCGTAGAAGTCACCGTTTGACTATCATGAATTCGCCAAAAGATCAAAGTAACTAAAGTATCACCAATGACGTGCC"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(1F, 1F, 2, 1), 1, 1, 6, 1, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "TCAGTGCCAAACCTACTCCTCAGGCAAATATTGAGCTCTTCTAGGAACTCTTTAACAGCAGTAGTGCAATCACCAGATATACACTCTGCAATGTAGAACGACACTTTTTCTAAGGTAGTTGTGAATCGTAAAGGTCACCTTTTCACTATCTTGAATTCACCAAAAGATCAAAGTAACCAAAGTATCACCAATGACGTGCC\r\n||||||.|||||||||..|||||||||||||.|||||.|.|||||||||||.|||||||||.|.|||||||||||||||||||||||||.|||||||||||||||||||.|||.||||.|||||||||||.|.||||||.|||.||||||.|||||||.||||||||||||||||||.||||||||||||||||||||||\r\nTCAGTGACAAACCTACGGCTCAGGCAAATATAGAGCTTTACTAGGAACTCTATAACAGCAGCATTGCAATCACCAGATATACACTCTGCTATGTAGAACGACACTTTTTTTAAAGTAGCTGTGAATCGTAGAAGTCACCGTTTGACTATCATGAATTCGCCAAAAGATCAAAGTAACTAAAGTATCACCAATGACGTGCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 160);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 200);
            Assert.AreEqual(sml.End_Query, 200);
            Assert.AreEqual(sml.UsedDynamicBanding, true);


            database = new Sequence_ByteArray("AGTGGATTTCAGCTAACCCACTGTCTTTTTTTTGGATACTACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGCAGTAGGGTATACATCTGATTACATTAAAGTGCGAATACTCTTTTGGTGGAGTTACCTCTAATTCCTTGTAGGGGGTTGGTATATGTAATACACGAGTTCAACTCTTCTTGGATCC"); //Start 14
            query = new Sequence_ByteArray("AGTGTATTTCCGCAAACCCACTGTCTTTGTTTTAAGATACCACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGTAGTTGGGTATACATCTGATTACATTAAAGTGTCAATACGCTTGTGATGGAGCTACCCCTGATCTTGTAGGGGGCTAGGACATGTCGTACACGAGTTCAACTCTTCTTTGATCC"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(), 1, 1, 20, 5, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AGTGGATTTCAGCTAACCCACTGTCTTTTTTTT-GGATACTACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGCAGTAGGGTATACATCTGATTACATTAAAGTGCGAATACTCTTTTGGTGGAGTTACCTCTAATTCCTTGTAGGGGGTTGGTATATGTAATACACGAGTTCAACTCTTCTTGGATCC\r\n||||.|||||.||.||||||||||||||.|||| .|||||.||||||||||||||||||||||||||||||||||||||||||||.|||.|||||||||||||||||||||||||||..|||||.|||.||.|||||.||||.||.||  |||||||||||.|.|.|.||||..|||||||||||||||||||||.|||||\r\nAGTGTATTTCCGCAAACCCACTGTCTTTGTTTTAAGATACCACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGTAGTTGGGTATACATCTGATTACATTAAAGTGTCAATACGCTTGTGATGGAGCTACCCCTGAT--CTTGTAGGGGGCTAGGACATGTCGTACACGAGTTCAACTCTTCTTTGATCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 227.7199) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 200);
            Assert.AreEqual(sml.End_Query, 199);
            Assert.AreEqual(sml.UsedDynamicBanding, true);


            database = new Sequence_ByteArray("AGTGGATTTCAGCTAACCCACTGTCTTTTTTTTGGATACTACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGCAGTAGGGTATACATCTGATTACATTAAAGTGCGAATACTCTTTTGGTGGAGTTACCTCTAATTCCTTGTAGGGGGTTGGTATATGTAATACACGAGTTCAACTCTTCTTGGATCC"); //Start 14
            query = new Sequence_ByteArray("AGTGTATTTCCGCAAACCCACTGTCTTTGTTTTAAGATACCACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGTAGTTGGGTATACATCTGATTACATTAAAGTGTCAATACGCTTGTGATGGAGCTACCCCTGATCTTGTAGGGGGCTAGGACATGTCGTACACGAGTTCAACTCTTCTTTGATCC"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(), 1, 1, 50, 20, 0);
            al = sml.GetAlignment();
            al.AlignmentFormater = new PairwiseAlignmentFormater_Blast();
            Assert.AreEqual(al.ToString(), "AGTGGATTTCAGCTAACCCACTGTCTTT-TTTTTGGATACTACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGCAGTAGGGTATACATCTGATTACATTAAAGTGCGAATACTCTTTTGGTGGAGTTACCTCTAATTCCTTGTAGGGGGTTGGTATATGTAATACACGAGTTCAACTCTTCTTGGATCC\r\n|||| ||||| || |||||||||||||| ||||  ||||| |||||||||||||||||||||||||||||||||||||||||||| ||| |||||||||||||||||||||||||||  ||||| ||| || ||||| |||| || ||  ||||||||||| | | | ||||  ||||||||||||||||||||| |||||\r\nAGTGTATTTCCGCAAACCCACTGTCTTTGTTTTAAGATACCACCAGCTGCGTTATTCGAGTATGAGGATTTTCCAACGTATCAAGTAGTTGGGTATACATCTGATTACATTAAAGTGTCAATACGCTTGTGATGGAGCTACCCCTGAT--CTTGTAGGGGGCTAGGACATGTCGTACACGAGTTCAACTCTTCTTTGATCC\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 227.7199) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 200);
            Assert.AreEqual(sml.End_Query, 199);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            database = new Sequence_ByteArray("AAATTCCCGGTTTAAGGGCCAAATT"); //Start 14
            query = new Sequence_ByteArray("AAAATCCCCGTTTTAGGGGCAAAAT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(5, 1), 1, 1, 5, 2, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAATTCCCGGTTTAAGGGCCAAA\r\n|||.||||.||||.||||.||||\r\nAAAATCCCCGTTTTAGGGGCAAA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 21.7F);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 23);
            Assert.AreEqual(sml.End_Query, 23);
            Assert.AreEqual(sml.UsedDynamicBanding, true);



            database = new Sequence_ByteArray("AAATTCCCGGTTTAAGGGCCAAATT"); //Start 14
            query = new Sequence_ByteArray("AAAATCCCCGTTTTAGGGGCAAAAT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(5, 1), 1, 1, 10, 5, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAATTCCCGGTTTAAGGGCCAAA\r\n|||.||||.||||.||||.||||\r\nAAAATCCCCGTTTTAGGGGCAAA\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 21.7F) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 23);
            Assert.AreEqual(sml.End_Query, 23);
            Assert.AreEqual(sml.UsedDynamicBanding, true);


            database = new Sequence_ByteArray("AAATTCCCGGTTTAAGGGCCAAATT"); //Start 14
            query = new Sequence_ByteArray("AAAATCCCCGTTTTAGGGGCAAAAT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(5, 1), 25, 25, 3, 1, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "AAATTCCCGGTTTAAGGGCCAAATT\r\n|||.||||.||||.||||.||||.|\r\nAAAATCCCCGTTTTAGGGGCAAAAT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 20.9F) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 25);
            Assert.AreEqual(sml.End_Query, 25);
            Assert.AreEqual(sml.UsedDynamicBanding, true);


            database = new Sequence_ByteArray("ACGTACGTACGTACGTACGTACGTACGT"); //Start 14
            query = new Sequence_ByteArray("ACGACGTACGACGTACGACGTACGT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(3, 1), 1, 1, 7, 3, 0);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGTACGTACGTACGTACGT\r\n||| ||||||| ||||||| ||||||||\r\nACG-ACGTACG-ACGTACG-ACGTACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 32.5F) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 28);
            Assert.AreEqual(sml.End_Query, 25);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            database = new Sequence_ByteArray("ACGTACGTACGTACGTACGTACGTACGT"); //Start 14
            query = new Sequence_ByteArray("ACGACGTACGACGTACGACGTACGT"); //Start 4
            sml = new AnchoredDynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix_PAM25WithDifferentlyScoredTransitionTransversion(3, 1), 1, 1, 7, 3, 100);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "ACGTACGTACGTACGTACGTACGTACGT\r\n||| ||||||| ||||||| ||||||||\r\nACG-ACGTACG-ACGTACG-ACGTACGT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.IsTrue(Math.Abs(sml.Score - 32.5F) < 0.001);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 28);
            Assert.AreEqual(sml.End_Query, 25);
            Assert.AreEqual(sml.UsedDynamicBanding, false);
        }

        [Test]
        public void Test_AnchoredDynamicBanded454SmithWatermanGotoh()
        {
            ISequenceContainer database = new Sequence_ByteArray("T");
            ISequenceContainer query = new Sequence_ByteArray("T");
            Anchored454DynamicBandedSmithWatermanGotoh sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), true, 1, 1, 3, 1, 0, 3);
            PairwiseAlignment al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 3.0);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), false, 1, 1, 3, 1, 0, 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 3.0);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

             database = new Sequence_ByteArray("T");
             query = new Sequence_ByteArray("T");
             sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), true, 1, 1, 3, 1, 2, 3);
             al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 3.0);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.UsedDynamicBanding, false);


            sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), false, 1, 1, 3, 1, 2, 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), "T\r\n|\r\nT\r\n");
            Assert.AreEqual(sml.Start_Database, 1);
            Assert.AreEqual(sml.Score, 3.0);
            Assert.AreEqual(sml.Start_Query, 1);
            Assert.AreEqual(sml.End_Database, 1);
            Assert.AreEqual(sml.End_Query, 1);
            Assert.AreEqual(sml.UsedDynamicBanding, false);


            //
            //Test nix
            //
            database = new Sequence_ByteArray("C");
            query = new Sequence_ByteArray("T");
            sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), false, 1, 1, 3, 1, 0, 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al,null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Score, 0.0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), true, 1, 1, 3, 1, 0, 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al, null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Score, 0.0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.UsedDynamicBanding, true);

            ///
            //Test direction
            //
            database = new Sequence_ByteArray("CGTCGTAAA");
            query = new Sequence_ByteArray("CAAAACAACAGTAAACAAAACAAAAACACA");
            sml = new Anchored454DynamicBandedSmithWatermanGotoh(database, query, SubstitutionMatrixFactory.GetNucleotideSequenceSubstitutionMatrix(3, 5, 11, 2.0F), true, 1, 1, 7, 2, 0, 3);
            al = sml.GetAlignment();
            Assert.AreEqual(al.ToString(), null);
            Assert.AreEqual(sml.Start_Database, 0);
            Assert.AreEqual(sml.Score, 0.0);
            Assert.AreEqual(sml.Start_Query, 0);
            Assert.AreEqual(sml.End_Database, 0);
            Assert.AreEqual(sml.End_Query, 0);
            Assert.AreEqual(sml.UsedDynamicBanding, true);


        }

    }



}


