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



namespace TestBioSeq
{
    using Bio.Seq;
    using Bio.Seq.IO;
    using Bio.Seq.Misc;
    using Bio.Seq.Sort;
    using Bio.Seq.Statistics;
    using Bio.Stat;

    using NUnit.Framework;
    

    public class Program
    {
        public static void Main()
        {
            Seq_Basic.Start();
            //Seq_IO.Start();
            Seq_Misc.Start();
            Seq_Sort.Start();
            Seq_Statistics.Start();
        }
    }

    /// <summary>
    /// Test the basic classes of Bio Seq
    /// </summary>
    [TestFixture]
    public class Seq_Basic
    {

        public static void Start()
        {
            Seq_Basic p = new Seq_Basic();
            p.Test_DNASequenceValidator_ATCG();
            p.Test_DNASequenceValidator_ATCGIndel();
            p.Test_DNASequenceValidator_ATCGN();
            p.Test_DNASequenceValidator_ATCGNYRWSKMBDHV();
            p.Test_DNASequenceValidator_UppercaseATCGNYRWSKMBDHV();
            p.Test_RNASequenceValidator_ACGU();
            p.Test_RNASequenceValidator_ACGUN();
            p.Test_Sequence_ByteArray();
            p.Test_Sequence_File();
            p.Test_Sequence_MergedString();
            p.Test_Sequence_String();
            p.Test_NucleotideSequence();
            p.Test_NucSeqInfo();
            p.Test_SequenceUtility();
            p.Test_NucleotideSequenceUtility();
            p.Test_DNADictionary();
            p.Test_HashATCG();
        }

        [Test]
        public void Test_DNASequenceValidator_ATCG()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val=new DNASequenceValidator_ATCG();

            string_valid= "AaAaaTTttTCCccCGggGGaaaaaTTTTTgggggCCCCC";
            string_invalid = "AAAAATTTTTCCCCC-----NNNNNYYYYYXXXXXGGGGGUUUUU_____";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach(char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 16);
            Assert.AreEqual(invalidPos.Count, 30);

 
        }

        [Test]
        public void Test_DNASequenceValidator_ATCGIndel()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new DNASequenceValidator_ATCGIndel();

            string_valid = "AaAaaTTttTCCccCGggGG-----aaaaaTTTTTgggggCCCCC";
            string_invalid = "AAAAATTTTTCCCCC-----NNNNNYYYYYXXXXXGGGGGUUUUU";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 21);
            Assert.AreEqual(invalidPos.Count, 20);


        }

        [Test]
        public void Test_DNASequenceValidator_ATCGN()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new DNASequenceValidator_ATCGN();

            string_valid = "AaAaaTTttTCCccCGggGGnnnNNNaaaaaTTTTTgggggCCCCC";
            string_invalid = "AAAAATTTTTCCCCC-----NNNNNYYYYYXXXXXGGGGGUUUUU";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 16);
            Assert.AreEqual(invalidPos.Count, 20);
        }

        [Test]
        public void Test_DNASequenceValidator_ATCGNYRWSKMBDHV()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new DNASequenceValidator_ATCGNYRWSKMBDHV();

            string_valid = "AaAaaTTttTCCccCGggGGnnnNNNaaaaaTTTTTgggggCCCCCYYyyYRRRRRWWwwwssssskkkkkmmmmmbbbbbdddddhhhhhvvvvv";
            string_invalid = "AAAAATTTTTgggggCCCCC-----NNNNNYYyyYRRRRRWWwwwssssskkkkkmmmmmbbbbbdddddhhhhhvvvvv";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 21);
            Assert.AreEqual(invalidPos.Count, 5);
        }

        [Test]
        public void Test_DNASequenceValidator_UppercaseATCGNYRWSKMBDHV()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new DNASequenceValidator_UppercaseATCGNYRWSKMBDHV();

            string_valid = "AAAAATTTTTCCCCCGGGGGNNNNNAAAAATTTTTCCCCCYYYYYRRRRRWWWWWWSSSSSKKKKKMMMMMBBBBBBDDDDDHHHHHVVVVV";
            string_invalid = "aaaaaTTTTTCCCCCGGGGGNNNNN-----AAAAATTTTTCCCCCYYYYYUUUUURRRRRWWWWWWSSSSSKKKKKMMMMMBBBBBBDDDDDHHHHHVVVVV";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 1);
            Assert.AreEqual(invalidPos.Count, 15);
        }

        [Test]
        public void Test_RNASequenceValidator_ACGU()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new RNASequenceValidator_ACGU();

            string_valid = "AAAACCCCGGGGGCAUUUUU";
            string_invalid = "AAAAACCCCCGGGGGTTTTTUUUUU-----YYYYY";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 16);
            Assert.AreEqual(invalidPos.Count, 15);
        }

        [Test]
        public void Test_RNASequenceValidator_ACGUN()
        {
            ISequenceValidator val;
            string string_valid;
            Sequence_String sequence_valid;
            NucleotideSequence ns_valid;
            INucSeqInfo nsi_valid;

            string string_invalid;
            Sequence_String sequence_invalid;
            NucleotideSequence ns_invalid;
            INucSeqInfo nsi_invalid;
            List<int> invalidPos;

            val = new RNASequenceValidator_ACGUN();

            string_valid = "AAAACCCCGGGGGCAUUUUUNNNNN";
            string_invalid = "AAAAACCCCCGNNNNNGGGGTTTTTUUUUU-----YYYYY";
            sequence_valid = new Sequence_String(string_valid);
            sequence_invalid = new Sequence_String(string_invalid);
            ns_valid = new NucleotideSequence("hans", sequence_valid);
            ns_invalid = new NucleotideSequence("hans", sequence_invalid);
            nsi_valid = ns_valid.GetNucSeqInfo();
            nsi_invalid = ns_invalid.GetNucSeqInfo();

            Assert.IsTrue(val.IsValid(sequence_valid));
            Assert.IsFalse(val.IsValid(sequence_invalid));

            Assert.IsTrue(val.IsValid(nsi_valid));
            Assert.IsFalse(val.IsValid(nsi_invalid));

            foreach (char c in string_valid)
            {
                Assert.IsTrue(val.IsValid(c));
            }

            invalidPos = val.GetInvalid(sequence_valid);
            Assert.IsEmpty(invalidPos);
            invalidPos = val.GetInvalid(sequence_invalid);
            Assert.IsNotEmpty(invalidPos);
            Assert.AreEqual(invalidPos[0], 21);
            Assert.AreEqual(invalidPos.Count, 15);
        }

        [Test]
        public void Test_Sequence_ByteArray()
        {
            ISequenceContainer seqC;
            ISequenceContainer subC;
            string test = "";
            string nt = "ACGTNACGTN";
            char[] ntc = nt.ToCharArray();
            byte[] ntb = new byte[ntc.Length];
            for (int i = 0; i < ntc.Length; i++)
            {
                ntb[i] = (byte)ntc[i];
            }
            seqC = new Sequence_ByteArray("ACGTNACGTN");
            seqC.Append(nt);
            seqC.Append(ntc);
            seqC.Append(ntb);
            foreach (char c in nt)
            {
                seqC.Append((byte)c);
                seqC.Append(c);
            }



            Assert.AreEqual(60, seqC.Length);
            Assert.AreEqual(seqC.ToString(), "ACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNAACCGGTTNNAACCGGTTNN");
            Assert.IsTrue(new DNASequenceValidator_ATCGN().IsValid(seqC));
            Assert.IsTrue(seqC.NewSequence() is Sequence_ByteArray);
            Assert.IsTrue(seqC.ConvertTo(new Sequence_String()) is Sequence_String);


            for (int i = 0; i < 60; i++)
            {
                test += seqC[i];
            }
            Assert.AreEqual(seqC.ToString(), test);
            subC = seqC.SubSequence(10, 20);
            Assert.IsTrue(subC is Sequence_ByteArray);
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", subC.ToString());
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", seqC.ToString(10, 20));
            seqC.Clear();
            Assert.AreEqual(seqC.Length, 0);
            seqC.Append("aaaaatttttcccccgggggaaaaatttttcccccggggg");
            Assert.AreEqual(seqC.ToString(5, 10), "tttttccccc");
            seqC.ToUpper();
            Assert.AreEqual(seqC.ToString(5, 10), "TTTTTCCCCC");


            //Test TrimExcess
            Assert.AreEqual(seqC.Capacity, 64);
            seqC.TrimExcess();
            Assert.AreEqual(seqC.Capacity, 40);

        }

        [Test]
        public void Test_Sequence_File()
        {
            ISequenceContainer seqC;
            ISequenceContainer subC;
            string test = "";
            string nt = "ACGTNACGTN";
            string fullString="ACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNAACCGGTTNNAACCGGTTNN";
            char[] ntc = nt.ToCharArray();
            byte[] ntb = new byte[ntc.Length];
            for (int i = 0; i < ntc.Length; i++)
            {
                ntb[i] = (byte)ntc[i];
            }
            seqC = new Sequence_File("ACGTNACGTN");
            seqC.Append(nt);
            seqC.Append(ntc);
            seqC.Append(ntb);
            foreach (char c in nt)
            {
                seqC.Append((byte)c);
                seqC.Append(c);
            }
            Assert.AreEqual(60, seqC.Length);
            Assert.AreEqual(seqC.ToString(),fullString );
            for (int i = 0; i < fullString.Length; i++)
            {
                Assert.AreEqual(seqC[i], fullString[i]);
            }
            Assert.IsTrue(new DNASequenceValidator_ATCGN().IsValid(seqC));
            Assert.IsTrue(seqC.NewSequence() is Sequence_File);
            Assert.IsTrue(seqC.ConvertTo(new Sequence_String()) is Sequence_String);
            for (int i = 0; i < 60; i++)
            {
                test += seqC[i];
            }
            Assert.AreEqual(seqC.ToString(), test);
            subC = seqC.SubSequence(10, 20);
            Assert.IsTrue(subC is Sequence_File);
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", subC.ToString());
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", seqC.ToString(10, 20));

            //Clear
            seqC.Clear();
            Assert.AreEqual(seqC.Length, 0);
            seqC.Append("aaaaatttttcccccgggggaaaaatttttcccccggggg");
            Assert.AreEqual(seqC.ToString(5, 10), "tttttccccc");
            seqC.ToUpper();
            Assert.AreEqual(seqC.ToString(5, 10), "TTTTTCCCCC");
            //Test TrimExcess
            Assert.AreEqual(seqC.Capacity, 40);
            seqC.TrimExcess();
            Assert.AreEqual(seqC.Capacity, 40);
            seqC[1] = 'C';

            Assert.AreEqual(seqC.ToString(), "ACAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG");


        }

        [Test]
        public void Test_Sequence_MergedString()
        {
            ISequenceContainer seqC;
            ISequenceContainer subC;
            string test = "";
            string nt = "ACGTNACGTN";
            char[] ntc = nt.ToCharArray();
            byte[] ntb = new byte[ntc.Length];
            for (int i = 0; i < ntc.Length; i++)
            {
                ntb[i] = (byte)ntc[i];
            }
            seqC = new Sequence_MergedString("ACGTNACGTN");
            seqC.Append(nt);
            seqC.Append(ntc);
            seqC.Append(ntb);
            foreach (char c in nt)
            {
                seqC.Append((byte)c);
                seqC.Append(c);
            }



            Assert.AreEqual(60, seqC.Length);
            Assert.AreEqual(seqC.ToString(), "ACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNAACCGGTTNNAACCGGTTNN");
            Assert.IsTrue(new DNASequenceValidator_ATCGN().IsValid(seqC));
            Assert.IsTrue(seqC.NewSequence() is Sequence_MergedString);
            Assert.IsTrue(seqC.ConvertTo(new Sequence_ByteArray()) is Sequence_ByteArray);


            for (int i = 0; i < 60; i++)
            {
                test += seqC[i];
            }
            Assert.AreEqual(seqC.ToString(), test);
            subC = seqC.SubSequence(10, 20);
            Assert.IsTrue(subC is Sequence_MergedString);
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", subC.ToString());
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", seqC.ToString(10, 20));
            seqC.Clear();
            Assert.AreEqual(seqC.Length, 0);
            seqC.Append("aaaaatttttcccccgggggaaaaatttttcccccggggg");
            Assert.AreEqual(seqC.ToString(5, 10), "tttttccccc");
            seqC.ToUpper();
            Assert.AreEqual(seqC.ToString(5, 10), "TTTTTCCCCC");


        }

        [Test]
        public void Test_Sequence_String()
        {
            ISequenceContainer seqC;
            ISequenceContainer subC;
            string test = "";
            string nt = "ACGTNACGTN";
            char[] ntc = nt.ToCharArray();
            byte[] ntb = new byte[ntc.Length];
            for (int i = 0; i < ntc.Length; i++)
            {
                ntb[i] = (byte)ntc[i];
            }
            seqC = new Sequence_String("ACGTNACGTN");
            seqC.Append(nt);
            seqC.Append(ntc);
            seqC.Append(ntb);
            foreach (char c in nt)
            {
                seqC.Append((byte)c);
                seqC.Append(c);
            }



            Assert.AreEqual(60, seqC.Length);
            Assert.AreEqual(seqC.ToString(), "ACGTNACGTNACGTNACGTNACGTNACGTNACGTNACGTNAACCGGTTNNAACCGGTTNN");
            Assert.IsTrue(new DNASequenceValidator_ATCGN().IsValid(seqC));
            Assert.IsTrue(seqC.NewSequence() is Sequence_String);
            Assert.IsTrue(seqC.ConvertTo(new Sequence_ByteArray()) is Sequence_ByteArray);


            for (int i = 0; i < 60; i++)
            {
                test += seqC[i];
            }
            Assert.AreEqual(seqC.ToString(), test);
            subC = seqC.SubSequence(10, 20);
            Assert.IsTrue(subC is Sequence_String);
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", subC.ToString());
            Assert.AreEqual("ACGTNACGTNACGTNACGTN", seqC.ToString(10, 20));
            seqC.Clear();
            Assert.AreEqual(seqC.Length, 0);
            seqC.Append("aaaaatttttcccccgggggaaaaatttttcccccggggg");
            Assert.AreEqual(seqC.ToString(5, 10), "tttttccccc");
            seqC.ToUpper();
            Assert.AreEqual(seqC.ToString(5, 10), "TTTTTCCCCC");
        }

        [Test]
        public void Test_NucleotideSequence()
        {
            ISequenceContainer seq;
            string name;
            string parentName;
            string comment;
            string seqString;
                
                
            seqString = "AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNYYYYYYYYYYRRRRRRRRRRWWWWWWWWWWSSSSSSSSSS" +
    "KKKKKKKKKKMMMMMMMMMMBBBBBBBBBBDDDDDDDDDD----------HHHHHHHHHHVVVVVVVVVVUUUUUUUUUUXXXXXXXXXX";
            seq = new Sequence_ByteArray(seqString);
            name= "testsequence";
            parentName = "testsequenceparent";
            comment = "this is the test sequence which should be tested";

            NucleotideSequence test = new NucleotideSequence(name, parentName, comment, seq, null, null);

            Assert.AreEqual(test.Sequence.ToString(), seqString);
            Assert.AreEqual(test.Name, name);
            Assert.AreEqual(test.ParentName, parentName);
            Assert.AreEqual(test.Start, null);
            Assert.AreEqual(test.End, null);
            Assert.AreEqual(test.Comment, comment);
            Assert.AreEqual(test.Length, seqString.Length);
            Assert.AreEqual(test.ChildSequence, null);
            Assert.AreEqual(test.ParentSequence, null);

            Assert.AreEqual(test.Count_A, 10);
            Assert.AreEqual(test.Count_T, 10);
            Assert.AreEqual(test.Count_U, 10);
            Assert.AreEqual(test.Count_ATU, 30);
            Assert.AreEqual(test.Count_C, 10);
            Assert.AreEqual(test.Count_Indel, 10);
            Assert.AreEqual(test.Count_G, 10);
            Assert.AreEqual(test.Count_CG, 20);
            Assert.AreEqual(test.Count_Special, 100);
            Assert.AreEqual(test.Count_Invalid, 20);

            bool exceptionThrown = false;
            try
            {
                test.IsValid();
            }
            catch (Exception e)
            {
                Assert.IsTrue(e is ArgumentNullException);
                exceptionThrown = true;
            }
            Assert.IsTrue(exceptionThrown);


            test.Validator = new DNASequenceValidator_ATCGNYRWSKMBDHV();
            Assert.IsFalse(test.IsValid());
            test.Validator = new SequenceValidator_AllValid();
            Assert.IsTrue(test.IsValid());
            Assert.AreEqual(test.Length, 180);

            test = new NucleotideSequence(name, parentName, new Sequence_MergedString("AAAAAAAAAATTTTTTTTTTCCCCCGGGGGGGGGGCCCCC"), 20, 50);
            Assert.AreEqual(test.Start, 20);
            Assert.AreEqual(test.End, 50);
            test.IsRoot = false;
            Assert.AreEqual(test.IsRoot, false);
            Assert.IsTrue(Math.Abs(test.GC_Percent - 50.0) < 0.0001);

            //Test GC-content with invalid characters, invalid characters should just be ignored
            test = new NucleotideSequence(name, parentName, new Sequence_MergedString("AAAAAAAAAA----------XXXXXXXXXXTTTTTTTTTTCCCCCGGGGGGGGGGCCCCC"));
            Assert.IsTrue(Math.Abs(test.GC_Percent - 50.0) < 0.0001);


        }

        [Test]
        public void Test_NucSeqInfo()
        {
            ISequenceContainer seq;
            string name;
            string parentname;
            string comment;
            string seqString;

            seqString ="AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNYYYYYYYYYYRRRRRRRRRRWWWWWWWWWWSSSSSSSSSS" +
    "KKKKKKKKKKMMMMMMMMMMBBBBBBBBBBDDDDDDDDDD----------HHHHHHHHHHVVVVVVVVVVUUUUUUUUUUXXXXXXXXXX";


            seq = new Sequence_MergedString(seqString);
             name = "testsequenceaname";
             parentname = "testsequenceparentname";
             comment = "this is a comment about the testsequence string";
            int startPos = 50;
            int endPos = 100;
            NucleotideSequence nseq = new NucleotideSequence(name, parentname, comment, seq, startPos, endPos);
            nseq.Validator = new DNASequenceValidator_ATCGNYRWSKMBDHV();
            nseq.IsRoot = false;
            INucSeqInfo nsInfo = nseq.GetNucSeqInfo();

            Assert.AreEqual(nsInfo.Name, name);
            Assert.AreEqual(nsInfo.ParentName, parentname);
            Assert.AreEqual(nsInfo.Comment, comment);
            Assert.IsFalse(nsInfo.IsValid());
            Assert.AreEqual(nsInfo.IsRoot, false);
            Assert.AreEqual(nsInfo.Start, startPos);
            Assert.AreEqual(nsInfo.End, endPos);

            Assert.AreEqual(nsInfo.Length, seqString.Length);
            Assert.AreEqual(nsInfo.Count_A, 10);
            Assert.AreEqual(nsInfo.Count_T, 10);
            Assert.AreEqual(nsInfo.Count_U, 10);
            Assert.AreEqual(nsInfo.Count_ATU, 30);
            Assert.AreEqual(nsInfo.Count_C, 10);
            Assert.AreEqual(nsInfo.Count_G, 10);
            Assert.AreEqual(nsInfo.Count_CG, 20);
            Assert.AreEqual(nsInfo.Count_ATUCG, 50);
            Assert.AreEqual(nsInfo.Count_N, 10);
            Assert.AreEqual(nsInfo.Count_Special, 100);
            Assert.AreEqual(nsInfo.Count_Indel, 10);
            Assert.AreEqual(nsInfo.Count_Invalid_woIndel, 10);
            Assert.AreEqual(nsInfo.IsUpper, true);
            Assert.AreEqual(nsInfo.SequenceCount, 1);
            Assert.AreEqual(nsInfo.Length, nsInfo.Count_AllChar); //180

            Assert.AreNotSame(nsInfo, nsInfo.Copy());


            ///Round two
            NucSeqInfo ns2 = new NucSeqInfo("fritz", "blabla", 10, 20);
            ns2.Count_A = 10;
            ns2.Count_T = 2;
            ns2.Count_U = 3;
            ns2.Count_C = 4;
            ns2.Count_G = 5;
            ns2.Count_Indel = 10;
            ns2.Count_Invalid_woIndel = 6;
            ns2.Count_Special = 7; //37
            ns2.IsUpper = true;
            ns2.IsRoot = false;
            ns2.Comment = comment;
            ns2.Length = 47;
            ns2.ParentName = "franzferdinand";
            ns2.Validator = new DNASequenceValidator_ATCGNYRWSKMBDHV();

            NucSeqInfo sum = ns2 + ((NucSeqInfo)nsInfo);
            Assert.AreEqual(sum.ParentName, null);
            Assert.AreEqual(sum.Comment, comment);
            Assert.AreEqual(sum.Start, null);
            Assert.AreEqual(sum.End, null);
            Assert.AreEqual(sum.Count_A, 20);
            Assert.AreEqual(sum.Count_T, 12);
            Assert.AreEqual(sum.Count_C, 14);
            Assert.AreEqual(sum.Count_G, 15);
            Assert.AreEqual(sum.Count_Indel, 20);
            Assert.AreEqual(sum.Count_Special, 107);
            Assert.AreEqual(sum.Count_Invalid_woIndel, 16);
            Assert.AreEqual(sum.SequenceCount, 2);
            Assert.AreEqual(sum.Names.Count, 2);
            Assert.AreEqual(sum.Names[0], "fritz");
            Assert.AreEqual(sum.Names[1], "testsequenceaname");
            Assert.AreEqual(sum.IsUpper, true);
            Assert.AreEqual(sum.IsRoot, false);
            Assert.AreEqual(sum.Count_AllChar, 227);
            Assert.AreEqual(sum.Length, 227);
            NucSeqInfo f = new NucSeqInfo();
            NucSeqInfo s = new NucSeqInfo("robert");
            f = f + s;
            Assert.AreEqual(f.Names.Count, 1);
            Assert.AreEqual(f.Names[0], "robert");
            Assert.AreEqual(f.Name, "robert");



        }

        [Test]
        public void Test_SequenceUtility()
        {
            string dna = "aAtTcCgGnNyYxXaAtT";
            string dnaUpper = "AATTCCGGNNYYXXAATT";
            string rna = "aAuUcCgGnNyYxXaAuU";

            Assert.AreEqual(SequenceUtility.ConvertToRNA(new Sequence_String(dna)).ToString(), rna);
            Assert.AreEqual(SequenceUtility.ConvertToUpper(new Sequence_ByteArray(dna)).ToString(), dnaUpper);
            Assert.AreEqual(SequenceUtility.ConvertToDNA(new Sequence_ByteArray(rna)).ToString(), dna);
            Assert.IsFalse(SequenceUtility.IsUpper(new Sequence_String(rna)));
            Assert.IsTrue(SequenceUtility.IsUpper(new Sequence_String(dnaUpper)));

            //Translate into peptide chain
            string orf = "ATGUUUUUUATAAAACCCGGGCCUAGCUAGGUAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
            ISequenceContainer peptide = SequenceUtility.ConvertToPeptide(new Sequence_String(orf));


            //MASK INVALID chars
            ISequenceContainer masked = SequenceUtility.MaskSequence(new Sequence_ByteArray(dnaUpper), new DNASequenceValidator_ATCGN(), 'N');
            Assert.AreEqual(masked.ToString(), "AATTCCGGNNNNNNAATT");
            masked = new Sequence_ByteArray(dnaUpper);
            SequenceUtility.MaskThisSequence(masked, new DNASequenceValidator_ATCGN(), 'N');
            Assert.AreEqual(masked.ToString(), "AATTCCGGNNNNNNAATT");
            masked = SequenceUtility.MaskSequence(new Sequence_ByteArray(dnaUpper), new RNASequenceValidator_ACGUN(), 'N');
            Assert.AreEqual(masked.ToString(), "AANNCCGGNNNNNNAANN");

            ISequenceContainer seq = new Sequence_ByteArray("AAATTAATT"); //Length:9 stepSize:4 Hashcontent:6
            Dictionary<string, List<int>> hash = SequenceUtility.GetHashTable(seq, 4);
            Assert.AreEqual(hash.Count, 5);
            Assert.AreEqual(hash["AAAT"][0], 1);
            Assert.AreEqual(hash["AATT"].Count, 2);
            Assert.AreEqual(hash["AATT"][0], 2);
            Assert.AreEqual(hash["AATT"][1], 6);


        }

        [Test]
        public void Test_NucleotideSequenceUtility()
        {
            ISequenceContainer seq1 = new Sequence_ByteArray("ATCGNCCCCCCCCCCATATATATATATATATATATGGGGGGGGGGgggggggggg");
            string transSeq1 = "CCCCCCCCCCCCCCCCCCCCATATATATATATATATATATGGGGGGGGGGNCGAT";
            ISequenceContainer seq2 = new Sequence_ByteArray("YYYNNNBBBRRRWWWSSSKKKMMMDDDHHHVVVAATTCCGG");
            string transSeq2 = "CCGGAATTBBBDDDHHHKKKMMMSSSWWWYYYVVVNNNRRR";

            Assert.AreEqual(NucleotideSequenceUtility.GetReverseComplement(seq1).ToString(), transSeq1);
            Assert.AreEqual(NucleotideSequenceUtility.GetReverseComplement(seq2).ToString(), transSeq2);
            Assert.AreEqual(NucleotideSequenceUtility.GetReverseComplement(transSeq2), seq2.ToString());
            seq1.ToUpper();
            Assert.AreEqual(NucleotideSequenceUtility.GetReverseComplement(transSeq1), seq1.ToString());

        }

        [Test]
        public void Test_HashATCG()
        {
            Hash_ATCG h = new Hash_ATCG(1);
            Assert.AreEqual(h.GetHash("A"), 0);
            Assert.AreEqual(h.GetHash("T"), 1);
            Assert.AreEqual(h.GetHash("C"), 2);
            Assert.AreEqual(h.GetHash("G"), 3);
            Assert.AreEqual(h.GetHash("X"), null);

            h = new Hash_ATCG(2);
            Assert.AreEqual(h.GetHash("AA"), 0);
            Assert.AreEqual(h.GetHash("AT"), 1);
            Assert.AreEqual(h.GetHash("AC"), 2);
            Assert.AreEqual(h.GetHash("AG"), 3);
            Assert.AreEqual(h.GetHash("TA"), 4);
            Assert.AreEqual(h.GetHash("TT"), 5);
            Assert.AreEqual(h.GetHash("TC"), 6);
            Assert.AreEqual(h.GetHash("TG"), 7);
            Assert.AreEqual(h.GetHash("TN"), null);

            h = new Hash_ATCG(3);
            Assert.AreEqual(h.GetHash("TTT"), 21);

            h = new Hash_ATCG(4);
            Assert.AreEqual(h.GetHash("GGGG"), 255);
            Assert.AreEqual(h.GetHash("GGGC"), 254);
            Assert.AreEqual(h.GetHash("CAAA"), 128);
            Assert.AreEqual(h.GetHash("CANA"), null);

            h = new Hash_ATCG(8);
            Assert.AreEqual(h.GetHash("CAAAAAAA"), 32768);

            h = new Hash_ATCG(12);
            Assert.AreEqual(h.GetHash("CAAAAAAAAAAA"), 8388608);

            h = new Hash_ATCG(16);
            Assert.AreEqual(h.GetHash("CAAAAAAAAAAAAAAA"), 2147483648);

            h = new Hash_ATCG(20);
            Assert.AreEqual(h.GetHash("Caaaaaaaaaaaaaaaaaaa"), 549755813888);

            h = new Hash_ATCG(24);

            Assert.AreEqual(h.GetHash("caaaaaaaaaaaaaaaaaaaaaaa"), 140737488355328);
            Assert.AreEqual(h.GetHash("caaaaaaaaaaaaaaaaaaaaaat"), 140737488355329);
            Assert.AreEqual(h.GetHash("caaaaaaaaaaaaaaaaaaaaaac"), 140737488355330);
            Assert.AreEqual(h.GetHash("CCCCXXXXAAAATTTTTTTTAAAA"), null);
        }

        [Test]
        public void Test_DNADictionary()
        {
            ISequenceContainer s1 = new Sequence_ByteArray("AAT");
            ISequenceContainer s2 = new Sequence_ByteArray("ATA");
            ISequenceContainer s3 = new Sequence_ByteArray("ATN");
            DNADictionary<int> dict = new DNADictionary<int>(3);
            Assert.IsFalse(dict.ContainsKey(s3));
            Assert.IsTrue(dict.ContainsKey(s1));
            Assert.IsTrue(dict.ContainsKey(s2));
            dict[s1]++;
            dict[s2] += 3;
            Assert.AreEqual(dict[s1], 1);
            Assert.AreEqual(dict[s2], 3);
            Assert.AreEqual(dict[new Sequence_ByteArray("TTT")], 0);




        }

    }

    [TestFixture]
    public class Seq_IO
    {

        public static void Start()
        {
            Seq_IO p = new Seq_IO();
            p.Test_FastaReader_Sequence();
            p.Test_FastaWriter_Sequence();
            p.Test_FastaReader_NucSeqInfo();
            p.Test_FastaReaderChopper_NucleotideSequence();
            p.Test_ByteReader();
            p.Test_ByteWriter();
            p.Test_QualityFileReader();
            p.Test_QualityFileWriter();
            p.Test_NucleotideSequenceReaderDecorator_ObtainQualityScores();
            p.Test_NucSeqInfoReader();
            p.Test_NucSeqInfoWriter();
            p.Test_RandomGenerator_NucleotideSequence();
            p.Test_SequenceReaderDecorator_ExtractFlankingSequence();
            p.Test_SequenceReaderDecorater_UniqueSequenceID();
            p.Test_FastasReader_NucleotideSequence();

        }


        [Test]
        public void Test_FastaReader_Sequence()
        {
            string testFilePath = @"..\..\TestFiles\Nunit_seq_fastreader_nucleotideSequence.txt";
            string sequenceFirst = "AAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNNNNNNNNNNNNAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCC";
            string sequenceSecond = "UUUUUUUUUUUUUUUUUUUUUAAAAAAAAAAAAAAAAAAAAYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYUUUUUUUUUUNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            string commentSecond = "this is a comment";
            string sequenceThird = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNnnnnn";
            string sequenceThirdUppercase = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            string sequenceFourth = "123456789010111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182";
            string commentFourth = "again with comment";

            FastaReader_Sequence<NucleotideSequence> fr = new FastaReader_Sequence<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid());
            Assert.AreEqual(fr.SequencePrototype.GetType(), SequenceFactory.GetDefaultSequence().GetType());
            List<NucleotideSequence> nucs = new List<NucleotideSequence>();
            NucleotideSequence ns;
            int i = 0;
            while ((ns=fr.GetNextSequence())!=null)
            {
                nucs.Add(ns);
                i++;
            }
            fr.Close();
            Assert.AreEqual(i, 4);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[0].Sequence.ToString(), sequenceFirst);
            Assert.AreEqual(nucs[0].Length, 100);
            Assert.AreEqual(nucs[0].Sequence.Capacity, 128);
            Assert.AreEqual(nucs[1].Name, "second");
            Assert.AreEqual(nucs[1].Sequence.ToString(), sequenceSecond);
            Assert.AreEqual(nucs[1].Comment, commentSecond);
            Assert.AreEqual(nucs[2].Name, "third");
            Assert.AreEqual(nucs[2].Sequence.ToString(), sequenceThird);
            Assert.AreEqual(nucs[3].Name, "fourth");
            Assert.AreEqual(nucs[3].Sequence.ToString(), sequenceFourth);
            Assert.AreEqual(nucs[3].Comment, commentFourth);


            //Round TWO
            //Reads all sequences at once
            //Ignores invalid
            fr = new FastaReader_Sequence<NucleotideSequence>(testFilePath, new DNASequenceValidator_ATCGNYRWSKMBDHV(), FastaReaderOptions.IgnoreInvalid);
            nucs = fr.GetAllSequences();
            Assert.AreEqual(nucs.Count, 2);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[0].Sequence.ToString(), sequenceFirst);
            Assert.AreEqual(nucs[1].Name, "third");
            Assert.AreEqual(nucs[1].Sequence.ToString(), sequenceThird);
            fr.Close();

            //Round THREE
            //Throws exception upon invalid sequence?
            fr = new FastaReader_Sequence<NucleotideSequence>(testFilePath, new DNASequenceValidator_ATCGNYRWSKMBDHV());
            int seqCount=1;
            bool exception=false;
            try
            {
                while((ns=fr.GetNextSequence())!=null)
                {
                    seqCount++;
                }
            }
            catch(InvalidSequenceException e)
            {
                exception=true;
                Assert.IsTrue(e.SequenceValidator is DNASequenceValidator_ATCGNYRWSKMBDHV);
                Assert.AreEqual(e.InvalidPositions[0], 1);
            }
            fr.Close();
            Assert.IsTrue(exception);
            Assert.AreEqual(seqCount,2);


            //Round FOUR
            //Trims the sequence and reads the correct amount of sequence
            fr = new FastaReader_Sequence<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid(),FastaReaderOptions.TrimExcess|FastaReaderOptions.ConvertToUppercase);
            nucs.Clear();
            List<NucleotideSequence> temp;
            i = 0;
            while ((temp = fr.GetSequences(2)).Count > 0)
            {
                Assert.AreEqual(temp.Count, 2);
                nucs.AddRange(temp);
                i++;
                temp.Clear();


            }
            fr.Close();
            Assert.AreEqual(i, 2);
            Assert.AreEqual(nucs.Count, 4);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[0].Sequence.ToString(), sequenceFirst);
            Assert.AreEqual(nucs[0].Length, 100);
            Assert.AreEqual(nucs[0].Sequence.Capacity, 100);
            Assert.AreEqual(nucs[2].Sequence.ToString(), sequenceThirdUppercase);


            //Congratulation the heart, i.e the fasta reader works!!
        }

        [Test]
        public void Test_FastaWriter_Sequence()
        {
            string testFilePath = @"..\..\TestFiles\Nunit_seq_fastawriter_nucleotideSequence.txt";
            string sequenceFirst = "AAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNNNNNNNNNNNNAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCC";
            string sequenceSecond = "UUUUUUUUUUUUUUUUUUUUUAAAAAAAAAAAAAAAAAAAAYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYUUUUUUUUUUNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            string commentSecond = "this is a comment";
            string sequenceThird = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
            string sequenceFourth = "123456789010111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182";
            string commentFourth = "again with comment";
            System.IO.StreamWriter sw = new System.IO.StreamWriter(testFilePath, false);
            List<NucleotideSequence> sequences = new List<NucleotideSequence>();

            sequences.Add(new NucleotideSequence("first", new Sequence_String(sequenceFirst)));
            sequences.Add(new NucleotideSequence("second", commentSecond, new Sequence_String(sequenceSecond)));
            sequences.Add(new NucleotideSequence("third", new Sequence_String(sequenceThird)));
            FastaWriter<NucleotideSequence> faw = new FastaWriter<NucleotideSequence>(sw);
            faw.WriteSequences(sequences);
            faw.WriteSequence(new NucleotideSequence("fourth", commentFourth, new Sequence_String(sequenceFourth)));
            sw.Close();

            //Since the readers are already tested we use it to validate the writer
            System.IO.StreamReader sr = new System.IO.StreamReader(testFilePath);
            FastaReader_Sequence<NucleotideSequence> fr = new FastaReader_Sequence<NucleotideSequence>(sr, new SequenceValidator_AllValid());
            Assert.AreEqual(fr.SequencePrototype.GetType(), SequenceFactory.GetDefaultSequence().GetType());
            List<NucleotideSequence> nucs = fr.GetAllSequences();
            sr.Close();

            Assert.AreEqual(nucs.Count, 4);

            Assert.AreEqual(nucs[0].Sequence.ToString(), sequenceFirst);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[1].Sequence.ToString(), sequenceSecond);
            Assert.AreEqual(nucs[1].Name, "second");
            Assert.AreEqual(nucs[1].Comment, commentSecond);
            Assert.AreEqual(nucs[2].Sequence.ToString(), sequenceThird);
            Assert.AreEqual(nucs[2].Name, "third");
            Assert.AreEqual(nucs[3].Sequence.ToString(), sequenceFourth);
            Assert.AreEqual(nucs[3].Name, "fourth");
            Assert.AreEqual(nucs[3].Comment, commentFourth);
        }
   
        [Test]
        public void Test_FastaReader_NucSeqInfo()
        {
            string testFilePath = @"..\..\TestFiles\Nunit_seq_fastreader_nucseqInfo.txt";



            FastaReader_NucSeqInfo fr = new FastaReader_NucSeqInfo(testFilePath, new SequenceValidator_AllValid(), FastaReaderNucSeqInfoOptions.FullDetails);

            List<INucSeqInfo> nucs = new List<INucSeqInfo>();
            INucSeqInfo temp;

            int i = 0;
            while ((temp = fr.GetNextSequence()) !=null)
            {
                nucs.Add(temp);
                i++;
            }

            //DETAIL LEVEL: Full

            string commentFirst = "A20 T10 U10 C10 N50 G20 Invalid10 lower";
            Assert.AreEqual(i, 4);
            Assert.IsFalse(nucs[0].IsUpper.Value);
            Assert.AreEqual(nucs[0].Count_A, 20);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[0].Count_C, 10);
            Assert.AreEqual(nucs[0].Count_T, 10);
            Assert.AreEqual(nucs[0].Count_U, 10);
            Assert.AreEqual(nucs[0].Count_G, 20);

            Assert.AreEqual(nucs[0].Count_Invalid_woIndel, 10);
            Assert.AreEqual(nucs[0].Count_N, 50);
            Assert.AreEqual(nucs[0].Comment, commentFirst);
            Assert.AreEqual(nucs[0].Length, 140);
            Assert.AreEqual(nucs[0].Length, nucs[0].Count_AllChar);
            //Assert.IsTrue(nucs[0].Validator is SequenceValidator_AllValid);
            Assert.AreEqual(nucs[0].Start, null);
            Assert.AreEqual(nucs[0].End, null);
            Assert.AreEqual(nucs[0].ParentName, "");



            string commentSecond = "U30 A20 Special30 N70 upper";
            Assert.IsTrue(nucs[1].IsUpper.Value);
            Assert.AreEqual(nucs[1].Count_U, 30);
            Assert.AreEqual(nucs[1].Name, "second");
            Assert.AreEqual(nucs[1].Count_A, 20);
            Assert.AreEqual(nucs[1].Count_N, 70);
            Assert.AreEqual(nucs[1].Count_Special, 30);
            Assert.AreEqual(nucs[1].Length, 150);
            Assert.AreEqual(nucs[1].Comment, commentSecond);
            Assert.AreEqual(nucs[1].Start, null);
            Assert.AreEqual(nucs[1].End, null);
            Assert.AreEqual(nucs[1].ParentName, "");


            string commentThird = "N300 upper";
            Assert.IsTrue(nucs[2].IsUpper.Value);
            Assert.AreEqual(nucs[2].Name, "third");
            Assert.AreEqual(nucs[2].Count_N, 300);
            Assert.AreEqual(nucs[2].Length, 300);
            Assert.AreEqual(nucs[2].Start, null);
            Assert.AreEqual(nucs[2].End, null);
            Assert.AreEqual(nucs[2].ParentName, "");
            Assert.AreEqual(nucs[2].Comment, commentThird);
            nucs[2].Validator = new DNASequenceValidator_ATCGN();
            Assert.IsTrue(nucs[2].IsValid());


            string commentFourth = "Comment 1: Invalid50 A25 T25 Special50 lower";
            Assert.AreEqual(nucs[3].Comment, commentFourth);
            Assert.AreEqual(nucs[3].Name, "fourth");
            Assert.IsFalse(nucs[3].IsUpper.Value);
            Assert.AreEqual(nucs[3].Count_Invalid_woIndel, 50);
            Assert.AreEqual(nucs[3].Count_A, 25);
            Assert.AreEqual(nucs[3].Count_T, 25);
            Assert.AreEqual(nucs[3].Count_Special, 50);
            Assert.AreEqual(nucs[3].Length, 150);


            //DETAIL LEVEL: LOW

            fr = new FastaReader_NucSeqInfo(testFilePath);
            nucs.Clear();
            nucs = fr.GetAllSequences();


            Assert.AreEqual(nucs.Count, 4);
            Assert.AreEqual(nucs[0].Name, "first");
            Assert.AreEqual(nucs[0].Comment, commentFirst);
            Assert.AreEqual(nucs[0].Length, 140);
            Assert.AreEqual(nucs[0].IsUpper, null);
            Assert.AreEqual(nucs[0].Count_AllChar, 0);



            Assert.AreEqual(nucs[1].Name, "second");
            Assert.AreEqual(nucs[1].Comment, commentSecond);
            Assert.AreEqual(nucs[1].IsRoot, null);
            Assert.AreEqual(nucs[1].Count_A, 0);
            Assert.AreEqual(nucs[1].Count_T, 0);

            Assert.AreEqual(nucs[2].Name, "third");
            Assert.AreEqual(nucs[2].Comment, commentThird);
            Assert.AreEqual(nucs[2].IsUpper, null);
            nucs[2].Validator = new DNASequenceValidator_ATCGN();
            Assert.IsFalse(nucs[2].IsValid());
            Assert.AreEqual(nucs[2].Count_CG, 0);
            Assert.AreEqual(nucs[2].Count_C, 0);
            Assert.AreEqual(nucs[2].Count_G, 0);
            Assert.AreEqual(nucs[2].Count_Invalid_woIndel, 0);
            Assert.AreEqual(nucs[2].Count_Special, 0);


            Assert.AreEqual(nucs[3].Name, "fourth");
            Assert.AreEqual(nucs[3].Comment, commentFourth);
            Assert.AreEqual(nucs[3].Length, 150);


            //Test GetSequences(number) and ignoreInvalid
            nucs = new List<INucSeqInfo>();
            List<INucSeqInfo> tempList;

            fr = new FastaReader_NucSeqInfo(testFilePath,new DNASequenceValidator_ATCGNYRWSKMBDHV(),FastaReaderNucSeqInfoOptions.FullDetails|FastaReaderNucSeqInfoOptions.IgnoreInvalid);
            nucs.Clear();

            i = 0;
            while ((tempList = fr.GetSequences(1)).Count>0)
            {
                nucs.AddRange(tempList);
                tempList.Clear();
                i++;
            }

            Assert.AreEqual(i, 1);
            Assert.AreEqual(nucs.Count, 1);

            Assert.AreEqual(nucs[0].Name, "third");
            Assert.AreEqual(nucs[0].Comment, commentThird);
            Assert.AreEqual(nucs[0].IsUpper, true);
            Assert.AreEqual(nucs[0].Count_CG, 0);
            Assert.AreEqual(nucs[0].Count_C, 0);
            Assert.AreEqual(nucs[0].Count_G, 0);
            Assert.AreEqual(nucs[0].Count_N, 300);
            Assert.AreEqual(nucs[0].Count_Invalid_woIndel, 0);
            Assert.AreEqual(nucs[0].Count_Special, 0);


        }

        [Test]
        public void Test_FastaReaderChopper_NucleotideSequence()
        {
            string testFilePath = @"..\..\TestFiles\Nunit_seq_fastaReaderChopper_nucleotideSequence.txt";
            string seqFirst = "1234567890101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384";
            string commentFirst = "again with comment 160";


            FastaReaderSequenceChopper<NucleotideSequence> fac = new FastaReaderSequenceChopper<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid(), FastaReaderOptions.None, 160, 10);

            List<NucleotideSequence> seq = fac.GetAllSequences();

            Assert.AreEqual(seq[0].Name, "first");
            Assert.AreEqual(seq[0].Sequence.ToString(), seqFirst);
            Assert.AreEqual(seq[0].Comment, commentFirst);
            Assert.AreEqual(seq[0].Length, 160);
            Assert.AreEqual(fac.SequencePrototype.GetType(), SequenceFactory.GetDefaultSequence().GetType());



            //
            //Next round
            //
            fac = new FastaReaderSequenceChopper<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid(), FastaReaderOptions.None, 50, 10);
            fac.SequencePrototype = new Sequence_ByteArray();
            NucleotideSequence temp;
            seq.Clear();
            int i = 0;


            while ((temp = fac.GetNextSequence())!=null)
            {
                if (temp.Name == "first")
                {
                    seq.Add(temp);
                    i++;
                }

            }

            Assert.AreEqual(seq.Count, 4);
            Assert.AreEqual(i, 4);
            string chunkFirst = "12345678901011121314151617181920212223242526272829";
            string chunkSecond = "25262728293031323334353637383940414243444546474849";
            string chunkThird = "45464748495051525354555657585960616263646566676869";
            string chunkFourth = "6566676869707172737475767778798081828384";



            Assert.AreEqual(seq[0].Sequence.ToString(), chunkFirst);
            Assert.AreEqual(seq[0].Name, "first");
            Assert.AreEqual(seq[0].Comment, commentFirst);
            Assert.IsTrue(seq[0].Sequence is Sequence_ByteArray);

            Assert.AreEqual(seq[1].Sequence.ToString(), chunkSecond);
            Assert.AreEqual(seq[1].Name, "first");
            Assert.AreEqual(seq[1].Comment, commentFirst);

            Assert.AreEqual(seq[2].Sequence.ToString(), chunkThird);
            Assert.AreEqual(seq[2].Name, "first");
            Assert.AreEqual(seq[2].Comment, commentFirst);

            Assert.AreEqual(seq[3].Sequence.ToString(), chunkFourth);
            Assert.AreEqual(seq[3].Name, "first");
            Assert.AreEqual(seq[3].Comment, commentFirst);


            //Round three
            fac = new FastaReaderSequenceChopper<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid(), FastaReaderOptions.TrimExcess|FastaReaderOptions.IgnoreInvalid, 300, 10);
            seq.Clear();
            List<NucleotideSequence> tempList;

            string seqSecond = "1234567890101112131415161718192021222324252627282930313233343536373839404142434445464748495051525354555657585960616263646566676869707172737475767778798081828384858687888990919293949596979899100101021031041051061071081091101111121131141151161171181191";
            string commentSecond = "250";
            string seqThird = "123456789010111213141516171819202122232425262728293031323334353637383940414243444546474849505152535455565758596061626364656667686970717273747576777879808182838485868788899091929394959697989910010102103104105106107108109110111112113114115116117118119120121123124125126127128129123013013113213313413513";
            string commentThird = "300";
            i = 0;
            while ((tempList = fac.GetSequences(2)).Count > 0)
            {
                seq.AddRange(tempList);
                tempList.Clear();
                i++;
            }
            Assert.AreEqual(seq.Count, 3);
            Assert.AreEqual(i, 2);
            Assert.AreEqual(seq[0].Name, "first");
            Assert.AreEqual(seq[0].Sequence.ToString(), seqFirst);
            Assert.AreEqual(seq[0].Comment, commentFirst);

            Assert.AreEqual(seq[1].Name, "second");
            Assert.AreEqual(seq[1].Comment, commentSecond);
            Assert.AreEqual(seq[1].Sequence.ToString(), seqSecond);

            Assert.AreEqual(seq[2].Name, "third");
            Assert.AreEqual(seq[2].Comment, commentThird);
            Assert.AreEqual(seq[2].Sequence.ToString(), seqThird);

            //Round four

            fac = new FastaReaderSequenceChopper<NucleotideSequence>(testFilePath, new SequenceValidator_AllValid(), FastaReaderOptions.None, 50, 10);
            seq.Clear();
            i = 0;
            while ((tempList = fac.GetSequences(2)).Count > 0)
            {
                seq.AddRange(tempList);
                tempList.Clear();
                i++;
            }

            Assert.AreEqual(seq.Count, 15);
            Assert.AreEqual(i, 8);

            string thirdFirst = "12345678901011121314151617181920212223242526272829";
            string thirdSecond = "25262728293031323334353637383940414243444546474849";
            string thirdThird = "45464748495051525354555657585960616263646566676869";
            string thirdFourth = "65666768697071727374757677787980818283848586878889";
            string thirdFivth = "85868788899091929394959697989910010102103104105106";
            string thirdSixth = "3104105106107108109110111112113114115116117118119120121123124125126127128129123013013113213313413513";

            Assert.AreEqual(seq[9].Name, "third");
            Assert.AreEqual(seq[9].Comment, commentThird);
            Assert.AreEqual(seq[9].Sequence.ToString(), thirdFirst);

            Assert.AreEqual(seq[10].Name, "third");
            Assert.AreEqual(seq[10].Comment, commentThird);
            Assert.AreEqual(seq[10].Sequence.ToString(), thirdSecond);

            Assert.AreEqual(seq[11].Sequence.ToString(), thirdThird);
            Assert.AreEqual(seq[12].Sequence.ToString(), thirdFourth);
            Assert.AreEqual(seq[13].Sequence.ToString(), thirdFivth);

            Assert.AreEqual(seq[14].Sequence.ToString(), thirdSixth);
            Assert.AreEqual(seq[14].Name, "third");
            Assert.AreEqual(seq[14].Comment, commentThird);








        }

        [Test]
        public void Test_ByteReader()
        {

            string path = @"..\..\TestFiles\Nunit_IOSeq_ByteReader.fasta";
            ByteReader_Sequence<NucleotideSequence> br = new ByteReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCG(),FastaReaderOptions.IgnoreInvalid);
            
            NucleotideSequence temp;
            List<NucleotideSequence> ns = new List<NucleotideSequence>();
            int count = 0;
            while ((temp = br.GetNextSequence())!=null)
            {
                ns.Add(temp);
                count++;
            }

            Assert.AreEqual(count, 4);
            Assert.AreEqual(ns.Count, 4);
            Assert.AreEqual(ns[0].Name, "fritzi");
            Assert.AreEqual(ns[0].Comment, "first sequence");
            Assert.AreEqual(ns[0].Sequence.ToString(), "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG");
            Assert.AreEqual(ns[1].Name, "franz");
            Assert.AreEqual(ns[1].Comment, null);
            Assert.AreEqual(ns[1].Sequence.ToString(), "A");
            Assert.AreEqual(ns[3].Name, "franzFerdinand");
            Assert.AreEqual(ns[3].Comment, "my name is franz ferdinand");
            Assert.AreEqual(ns[3].Sequence.ToString(), "TCGT");


            br = new ByteReader_Sequence<NucleotideSequence>(path, new SequenceValidator_AllValid());

            ns = new List<NucleotideSequence>();
            List<NucleotideSequence> tempList;
            count = 0;

            while ((tempList = br.GetSequences(2)).Count > 0)
            {
                ns.AddRange(tempList);
                count++;
                tempList.Clear();
            }

            Assert.AreEqual(count, 3);
            Assert.AreEqual(ns.Count, 5);
            Assert.AreEqual(ns[0].Name, "fritzi");
            Assert.AreEqual(ns[0].Comment, "first sequence");
            Assert.AreEqual(ns[0].Sequence.ToString(), "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG");
            Assert.AreEqual(ns[1].Name, "franz");
            Assert.AreEqual(ns[1].Comment, null);
            Assert.AreEqual(ns[1].Sequence.ToString(), "A");
            Assert.AreEqual(ns[4].Name, "franzFerdinand");
            Assert.AreEqual(ns[4].Comment, "my name is franz ferdinand");
            Assert.AreEqual(ns[4].Sequence.ToString(), "TCGT");


            ///Round three
            ///

            br = new ByteReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCG(),FastaReaderOptions.TrimExcess|FastaReaderOptions.IgnoreInvalid|FastaReaderOptions.ConvertToUppercase);

            ns = new List<NucleotideSequence>();
            ns = br.GetAllSequences();
            Assert.AreEqual(ns.Count, 4);
            Assert.AreEqual(ns[0].Name, "fritzi");
            Assert.AreEqual(ns[0].Comment, "first sequence");
            Assert.AreEqual(ns[0].Sequence.ToString(), "AAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGGAAAAATTTTTCCCCCGGGGG");
            Assert.AreEqual(ns[1].Name, "franz");
            Assert.AreEqual(ns[1].Comment, null);
            Assert.AreEqual(ns[1].Sequence.ToString(), "A");
            Assert.AreEqual(ns[3].Name, "franzFerdinand");
            Assert.AreEqual(ns[3].Comment, "my name is franz ferdinand");
            Assert.AreEqual(ns[3].Sequence.ToString(), "TCGT");
            

            //TODO Test the ByteReader with repetitive sequence
            /* Todo Test ByteReader with repetive sequennces
            path = @"G:\Programs\00TestSequencen\Nunit_IOSeq_ByteReaderLogRep.fasta";
            ByteReader_Sequence<Bio.Repetitive.LogRepetitiveSequence> bir = new ByteReader_Sequence<Bio.Repetitive.LogRepetitiveSequence>(path, new SequenceValidator_AllValid(), false, 2);
            List<Bio.Repetitive.LogRepetitiveSequence> nis = bir.GetSequences();
            Assert.AreEqual(nis.Count, 1);
            Assert.AreEqual(nis[0].Name, "logg odd");
            Assert.AreEqual(nis[0].Comment, "each byte smaller than 254 is tested");
            ISequenceContainer s = nis[0].Sequence;

            for (int i = 0; i < 508; i++)
            {
                Assert.AreEqual(((byte)(i % 254)), (byte)s[i]);
            }
             * 
             */




        }

        [Test]
        public void Test_ByteWriter()
        {
            string path = @"..\..\TestFiles\Nunit_IOSeq_ByteWriter.fasta";

            System.IO.FileStream fs = new System.IO.FileStream(path, System.IO.FileMode.Create, System.IO.FileAccess.Write);
            ByteWriter<NucleotideSequence> bw = new ByteWriter<NucleotideSequence>(fs);

            List<NucleotideSequence> ns = new List<NucleotideSequence>();
            //Sequence one
            ISequenceContainer ba = new Sequence_ByteArray();
            for (int i = 0; i < 508; i++)
            {

                ba.Append((byte)(i % 254));
            }
            ns.Add(new NucleotideSequence("Test routine for ByteWriter", "each byte smaller than 254 is tested twice", ba));

            Sequence_ByteArray bar = new Sequence_ByteArray();
            for (int i = 0; i < 254; i++)
            {

                bar.Append((byte)(i % 254));
            }
            ns.Add(new NucleotideSequence("Test", null, bar));
            ns.Add(new NucleotideSequence("test2", "hallo", new Sequence_ByteArray("AAAATTTTAAAATTTTCCCCCGGGGG")));
            ns.Add(new NucleotideSequence("morgen morgen", null, new Sequence_ByteArray("TTTTTTTTTTTTTTTTTVVVVVVVVVVV")));

            //Write sequences and reset everything
            bw.WriteSequences(ns);
            fs.Close();
            ns = null;
            fs = null;
            bw = null;
            ba = null;




            ByteReader_Sequence<NucleotideSequence> br = new ByteReader_Sequence<NucleotideSequence>(path, new SequenceValidator_AllValid());
            ns = br.GetAllSequences();
            Assert.AreEqual(ns.Count, 4);
            Assert.AreEqual(ns[0].Name, "Test routine for ByteWriter");
            Assert.AreEqual(ns[0].Comment, "each byte smaller than 254 is tested twice");
            ISequenceContainer s = ns[0].Sequence;

            for (int i = 0; i < 508; i++)
            {
                Assert.AreEqual(((byte)(i % 254)), (byte)s[i]);
            }

            Assert.AreEqual(ns[1].Name, "Test");
            Assert.AreEqual(ns[1].Comment, null);
            s = ns[1].Sequence;

            for (int i = 0; i < 254; i++)
            {
                Assert.AreEqual(((byte)(i % 254)), (byte)s[i]);
            }

            Assert.AreEqual(ns[2].Name, "test2");
            Assert.AreEqual(ns[2].Comment, "hallo");
            Assert.AreEqual(ns[2].Sequence.ToString(), "AAAATTTTAAAATTTTCCCCCGGGGG");
            Assert.AreEqual(ns[3].Name, "morgen morgen");
            Assert.AreEqual(ns[3].Comment, null);
            Assert.AreEqual(ns[3].Sequence.ToString(), "TTTTTTTTTTTTTTTTTVVVVVVVVVVV");
        }

        [Test]
        public void Test_QualityFileReader()
        {
            string inputFile = @"..\..\TestFiles\Nunit_IOSeq_QualityFileReader.qul";
            System.IO.StreamReader sr = new System.IO.StreamReader(inputFile);
            QualityFileReader fr = new QualityFileReader(sr, new SequenceValidator_AllValid());

            int count = 0;
            QualitySequence temp;
            List<QualitySequence> qualFiles = new List<QualitySequence>();
            while ((temp = fr.GetNextSequence())!=null)
            {
                qualFiles.Add(temp);

                count++;
            }
            sr.Close();

            Assert.AreEqual(qualFiles.Count, 6);
            Assert.AreEqual(count, 6);
            Assert.AreEqual(qualFiles[0].Name, "fritz1");
            Assert.AreEqual(qualFiles[0].Sequence.ToString(), "\"!#\" \"#\"%!!\"!\"\"\"\"&\"\t\"");

            Assert.AreEqual(qualFiles[1].Name, "fritz2");
            Assert.AreEqual(qualFiles[3].Sequence.ToString(), " \f\r!\a\n\"!!\" \b!        \"  ! \" ");

            Assert.AreEqual(qualFiles[4].Name, "fritz5");
            Assert.AreEqual(qualFiles[5].Name, "fritz6");
            Assert.AreEqual(qualFiles[5].Sequence.ToString(), "\" \" \a\b   \" \r  \" \v!");

            //
            //Round 2
            //
            sr = new System.IO.StreamReader(inputFile);
            fr = new QualityFileReader(sr, new SequenceValidator_AllValid(),QualityFileReaderOptions.TrimExcess|QualityFileReaderOptions.IgnoreInvalid);
            qualFiles = fr.GetAllSequences();
            sr.Close();

            Assert.AreEqual(qualFiles[0].Name, "fritz1");
            Assert.AreEqual(qualFiles[1].Name, "fritz2");
            Assert.AreEqual(qualFiles[2].Name, "fritz3");
            Assert.AreEqual(qualFiles[3].Name, "fritz4");
            Assert.AreEqual(qualFiles[4].Name, "fritz5");
            Assert.AreEqual(qualFiles[5].Name, "fritz6");

            //
            //Round 3
            //
            qualFiles.Clear();
            sr = new System.IO.StreamReader(inputFile);
            fr = new QualityFileReader(sr, new SequenceValidator_AllValid());
            List<QualitySequence> tempList;
            count = 0;
            while ((tempList = fr.GetSequences(2)).Count > 0)
            {
                qualFiles.AddRange(tempList);
                tempList.Clear();
                count++;
            }
            sr.Close();

            Assert.AreEqual(count, 3);
            Assert.AreEqual(qualFiles[0].Name, "fritz1");
            Assert.AreEqual(qualFiles[1].Name, "fritz2");
            Assert.AreEqual(qualFiles[2].Name, "fritz3");
            Assert.AreEqual(qualFiles[3].Name, "fritz4");
            Assert.AreEqual(qualFiles[4].Name, "fritz5");
            Assert.AreEqual(qualFiles[5].Name, "fritz6");

        }
        
        [Test]
        public void Test_QualityFileWriter()
        {
            string inputFile = @"..\..\TestFiles\Nunit_IOSeq_QualityFileWriter.qul";
            System.IO.StreamWriter sw = new System.IO.StreamWriter(inputFile, false, Encoding.ASCII);
            QualityFileWriter fw = new QualityFileWriter(sw);
            List<QualitySequence> toWrite = new List<QualitySequence>();
            toWrite.Add(new QualitySequence("fritz1", new Sequence_ByteArray("\"!#\" \"#\"%!!\"!\"\"\"\"&\"\t\"")));
            toWrite.Add(new QualitySequence("fritz2", new Sequence_ByteArray(" \f\r!\a\n\"!!\" \b!        \"  ! \" ")));
            toWrite.Add(new QualitySequence("fritz3", new Sequence_ByteArray("\" \" \a\b   \" \r  \" \v!")));
            fw.WriteSequences(toWrite);
            sw.Close();

            List<QualitySequence> read;
            System.IO.StreamReader sr = new System.IO.StreamReader(inputFile);
            QualityFileReader qr = new QualityFileReader(sr, new SequenceValidator_AllValid());
            read = qr.GetAllSequences();

            Assert.AreEqual(read[0].Name, "fritz1");
            Assert.AreEqual(read[0].Sequence.ToString(), "\"!#\" \"#\"%!!\"!\"\"\"\"&\"\t\"");
            Assert.AreEqual(read[1].Name, "fritz2");
            Assert.AreEqual(read[1].Sequence.ToString(), " \f\r!\a\n\"!!\" \b!        \"  ! \" ");
            Assert.AreEqual(read[2].Name, "fritz3");
            Assert.AreEqual(read[2].Sequence.ToString(), "\" \" \a\b   \" \r  \" \v!");

            sr.Close();
        }
    
        [Test]
        public void Test_NucleotideSequenceReaderDecorator_ObtainQualityScores()
        {
            string pathFasta = @"..\..\TestFiles\NucleotideSequenceDecorator_ObtainQualityScores.fasta";
            string pathQual = @"..\..\TestFiles\NucleotideSequenceDecorator_ObtainQualityScores.qul";
            ISequenceReader<NucleotideSequence> ns = new FastaReader_Sequence<NucleotideSequence>(pathFasta, new DNASequenceValidator_ATCG());
            NucleotideSequenceReaderDecorator_ObtainQualityScores<NucleotideSequence> nqsr = new NucleotideSequenceReaderDecorator_ObtainQualityScores<NucleotideSequence>(ns, new string[] { pathQual },1,false);
           
            List<NucleotideSequence> temp;
            int count = 0;
            List<NucleotideSequence> nucs = new List<NucleotideSequence>();
            while ((temp = nqsr.GetSequences(1)).Count > 0)
            {
                nucs.AddRange(temp);
                count++;
                temp.Clear();
            }

            Assert.AreEqual(count, 3);
            Assert.AreEqual(nucs.Count, 3);
            Assert.AreEqual(nucs[0].Name, "hans");
            Assert.AreEqual(nucs[0].QualitySequence.Name, "hans");
            Assert.AreEqual(nucs[0].Sequence.Length, 60);
            Assert.AreEqual(nucs[0].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[0].Sequence[0], 'A');
            Assert.AreEqual(nucs[0].QualitySequence.Sequence[0], (char)1);
            Assert.AreEqual(nucs[0].Sequence[59], 'T');
            Assert.AreEqual(nucs[0].QualitySequence.Sequence[59], (char)60);

            Assert.AreEqual(nucs[1].Name, "sepp");
            Assert.AreEqual(nucs[1].QualitySequence.Name, "sepp");
            Assert.AreEqual(nucs[1].Sequence.Length, 60);
            Assert.AreEqual(nucs[1].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[1].Sequence[0], 'C');
            Assert.AreEqual(nucs[1].QualitySequence.Sequence[0], (char)6);
            Assert.AreEqual(nucs[1].Sequence[59], 'G');
            Assert.AreEqual(nucs[1].QualitySequence.Sequence[59], (char)65);

            Assert.AreEqual(nucs[2].Name, "robert");
            Assert.AreEqual(nucs[2].QualitySequence.Name, "robert");
            Assert.AreEqual(nucs[2].Sequence.Length, 60);
            Assert.AreEqual(nucs[2].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[2].Sequence[0], 'A');
            Assert.AreEqual(nucs[2].QualitySequence.Sequence[0], (char)11);
            Assert.AreEqual(nucs[2].Sequence[59], 'C');
            Assert.AreEqual(nucs[2].QualitySequence.Sequence[59], (char)70);


            //
            //Round 2
            //           
            nucs = null;
            ns = new FastaReader_Sequence<NucleotideSequence>(pathFasta, new DNASequenceValidator_ATCG());
            nqsr = new NucleotideSequenceReaderDecorator_ObtainQualityScores<NucleotideSequence>(ns, new string[] { pathQual }, 1, false);
            nucs = nqsr.GetAllSequences();

            Assert.AreEqual(nucs.Count, 3);
            Assert.AreEqual(nucs[0].Name, "hans");
            Assert.AreEqual(nucs[0].QualitySequence.Name, "hans");
            Assert.AreEqual(nucs[0].Sequence.Length, 60);
            Assert.AreEqual(nucs[0].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[0].Sequence[0], 'A');
            Assert.AreEqual(nucs[0].QualitySequence.Sequence[0], (char)1);
            Assert.AreEqual(nucs[0].Sequence[59], 'T');
            Assert.AreEqual(nucs[0].QualitySequence.Sequence[59], (char)60);

            Assert.AreEqual(nucs[1].Name, "sepp");
            Assert.AreEqual(nucs[1].QualitySequence.Name, "sepp");
            Assert.AreEqual(nucs[1].Sequence.Length, 60);
            Assert.AreEqual(nucs[1].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[1].Sequence[0], 'C');
            Assert.AreEqual(nucs[1].QualitySequence.Sequence[0], (char)6);
            Assert.AreEqual(nucs[1].Sequence[59], 'G');
            Assert.AreEqual(nucs[1].QualitySequence.Sequence[59], (char)65);

            Assert.AreEqual(nucs[2].Name, "robert");
            Assert.AreEqual(nucs[2].QualitySequence.Name, "robert");
            Assert.AreEqual(nucs[2].Sequence.Length, 60);
            Assert.AreEqual(nucs[2].QualitySequence.Sequence.Length, 60);
            Assert.AreEqual(nucs[2].Sequence[0], 'A');
            Assert.AreEqual(nucs[2].QualitySequence.Sequence[0], (char)11);
            Assert.AreEqual(nucs[2].Sequence[59], 'C');
            Assert.AreEqual(nucs[2].QualitySequence.Sequence[59], (char)70);





        }

        [Test]
        public void Test_NucSeqInfoReader()
        {
            string testFilePath = @"..\..\TestFiles\Nunit_Seq_NucSeqInfoReader.txt";
            System.IO.StreamReader sr = new System.IO.StreamReader(testFilePath, Encoding.ASCII);
            INucSeqInfoReader nr = new NucSeqInfoReader(sr);
            INucSeqInfo temp;
            List<INucSeqInfo> nucs = new List<INucSeqInfo>();
            int countReads = 0;
            while ((temp = nr.GetNextSequence())!=null)
            {
                nucs.Add(temp);
                countReads++;
            }
            sr.Close();
            Assert.AreEqual(countReads, 3);
            Assert.AreEqual(nucs.Count, 3);
            Assert.AreEqual(nucs[0].Names.Count, 4);
            Assert.AreEqual(nucs[0].Names[0], "testsequenceaname");
            Assert.AreEqual(nucs[0].Names[1], "franz");
            Assert.AreEqual(nucs[0].Names[2], "fritz");
            Assert.AreEqual(nucs[0].Names[3], "robert");
            Assert.AreEqual(nucs[0].Comment, "lallelu nur der mann im mond schaut zu");
            Assert.AreEqual(nucs[0].Length, 170);
            Assert.AreEqual(nucs[0].Count_A, 11);
            Assert.AreEqual(nucs[0].Count_T, 12);
            Assert.AreEqual(nucs[0].Count_U, 13);
            Assert.AreEqual(nucs[0].Count_C, 14);
            Assert.AreEqual(nucs[0].Count_G, 15);
            Assert.AreEqual(nucs[0].Count_N, 16);
            Assert.AreEqual(nucs[0].Count_Special, 107);
            Assert.AreEqual(nucs[0].Count_Invalid_woIndel, 18);
            Assert.AreEqual(nucs[0].Count_Indel, 10);
            Assert.AreEqual(nucs[0].SequenceCount, 19);
            Assert.AreEqual(nucs[0].ParentName, "kaiser franz der schlaue");
            Assert.AreEqual(nucs[0].Start, 50);
            Assert.AreEqual(nucs[0].End, 100);
            Assert.AreEqual(nucs[0].IsUpper, true);
            Assert.AreEqual(nucs[1].Comment, "");
            Assert.AreEqual(nucs[1].Names.Count, 4);
            Assert.AreEqual(nucs[1].SequenceCount, 10);
            Assert.AreEqual(nucs[1].Count_Indel, 10);
            Assert.AreEqual(nucs[2].Name, "robert der beste");
            Assert.AreEqual(nucs[2].Names.Count, 1);
            Assert.AreEqual(nucs[2].Comment, "lallelu nur der ro schaut zu");
            Assert.AreEqual(nucs[2].Count_Special, 300);
            Assert.AreEqual(nucs[2].Count_Indel, 10);


            ///
            ///Round 2
            ///
             sr = new System.IO.StreamReader(testFilePath, Encoding.ASCII);
             nr = new NucSeqInfoReader(sr);
             nucs = null;
             nucs = nr.GetAllSequences();
             Assert.AreEqual(nucs.Count, 3);
             Assert.AreEqual(nucs[0].Names.Count, 4);
             Assert.AreEqual(nucs[0].Names[0], "testsequenceaname");
             Assert.AreEqual(nucs[0].Names[1], "franz");
             Assert.AreEqual(nucs[0].Names[2], "fritz");
             Assert.AreEqual(nucs[0].Names[3], "robert");
             Assert.AreEqual(nucs[0].Comment, "lallelu nur der mann im mond schaut zu");
             Assert.AreEqual(nucs[0].Length, 170);
             Assert.AreEqual(nucs[0].Count_A, 11);
             Assert.AreEqual(nucs[0].Count_T, 12);
             Assert.AreEqual(nucs[0].Count_U, 13);
             Assert.AreEqual(nucs[0].Count_C, 14);
             Assert.AreEqual(nucs[0].Count_G, 15);
             Assert.AreEqual(nucs[0].Count_N, 16);
             Assert.AreEqual(nucs[0].Count_Special, 107);
             Assert.AreEqual(nucs[0].Count_Invalid_woIndel, 18);
             Assert.AreEqual(nucs[0].Count_Indel, 10);
             Assert.AreEqual(nucs[0].SequenceCount, 19);
             Assert.AreEqual(nucs[0].ParentName, "kaiser franz der schlaue");
             Assert.AreEqual(nucs[0].Start, 50);
             Assert.AreEqual(nucs[0].End, 100);
             Assert.AreEqual(nucs[0].IsUpper, true);
             Assert.AreEqual(nucs[1].Comment, "");
             Assert.AreEqual(nucs[1].Names.Count, 4);
             Assert.AreEqual(nucs[1].SequenceCount, 10);
             Assert.AreEqual(nucs[1].Count_Indel, 10);
             Assert.AreEqual(nucs[2].Name, "robert der beste");
             Assert.AreEqual(nucs[2].Names.Count, 1);
             Assert.AreEqual(nucs[2].Comment, "lallelu nur der ro schaut zu");
             Assert.AreEqual(nucs[2].Count_Special, 300);
             Assert.AreEqual(nucs[2].Count_Indel, 10);
        }

        [Test]
        public void Test_NucSeqInfoWriter()
        {

            //seq_ATUCGNYRWSKMBDHVX10
            ISequenceContainer seq = new Sequence_MergedString("AAAAAAAAAATTTTTTTTTTUUUUUUUUUUCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNYYYYYYYYYYRRRRRRRRRRWWWWWWWWWWSSSSSSSSSS"+
                "KKKKKKKKKKMMMMMMMMMMBBBBBBBBBBDDDDDDDDDDHHHHHHHHHHVVVVVVVVVV-----XXXXXXXXXX");
            string testFilePath = @"..\..\TestFiles\Nunit_Seq_NucSeqInfoWriter.txt";
            string name = "testsequenceaname";
            //string parentname = "testsequenceparentname";
            string comment = "this is a comment about the testsequence string";
            NucleotideSequence nseq = new NucleotideSequence(name, comment, seq);
            nseq.Validator = new DNASequenceValidator_ATCGNYRWSKMBDHV();
            nseq.IsRoot = false;
            INucSeqInfo nsInfo = nseq.GetNucSeqInfo();
            nsInfo.Names.Add("franz");
            nsInfo.SequenceCount = 10;
            nsInfo.Names.Add("fritz");
            nsInfo.Names.Add("robert");
            nsInfo.Comment = comment;
            nsInfo.Length = 170;
            nsInfo.ParentName = "kaiser franz der schlaue";
            System.IO.StreamWriter sw = new System.IO.StreamWriter(testFilePath, false, Encoding.ASCII);


            INucSeqInfoWriter swri = new NucSeqInfoWriter(sw);
            List<INucSeqInfo> nss = new List<INucSeqInfo>();
            nss.Add(nsInfo);
            nsInfo = new NucSeqInfo();
            nsInfo.Names.Add("hallo");
            nsInfo.Names.Add("hallo2");
            nsInfo.Count_Indel = 10;
            nsInfo.Count_G = 5;
            nsInfo.Count_C = 5;
            nsInfo.Count_A = 2;
            nsInfo.Comment = "this is a comment";
            nsInfo.Count_N = 7;
            nsInfo.Count_T = 5;
            nss.Add(nsInfo);
            Assert.AreEqual(nss.Count, 2);
            swri.WriteSeqInfos(nss);
            sw.Close();

            //
            //
            //
            System.IO.StreamReader sr = new System.IO.StreamReader(testFilePath);
            NucSeqInfoReader reader = new NucSeqInfoReader(sr);
            List<INucSeqInfo> temp;
            List<INucSeqInfo> nucSeq = new List<INucSeqInfo>();
            int count = 0;
            while ((temp = reader.GetSequences(1)).Count > 0)
            {
                nucSeq.AddRange(temp);
                count++;
                temp.Clear();
            }

            Assert.AreEqual(count, 2);
            Assert.AreEqual(nucSeq.Count, 2);
            Assert.AreEqual(nucSeq[0].Count_Indel, 5);
            Assert.AreEqual(nucSeq[0].Count_Invalid_woIndel, 10);
            Assert.AreEqual(nucSeq[0].Count_A, 10);
            Assert.AreEqual(nucSeq[0].Count_T, 10);
            Assert.AreEqual(nucSeq[0].Count_C, 10);
            Assert.AreEqual(nucSeq[0].Count_G, 10);
            Assert.Contains("franz", nucSeq[0].Names);
            Assert.Contains("robert", nucSeq[0].Names);
            Assert.Contains("testsequenceaname", nucSeq[0].Names);
            Assert.AreEqual(nucSeq[0].Comment, "this is a comment about the testsequence string");
            Assert.AreEqual(nucSeq[0].ParentName, "kaiser franz der schlaue");
            Assert.AreEqual(nucSeq[0].SequenceCount, 10);

            Assert.AreEqual(nucSeq[1].Count_Indel, 10);
            Assert.AreEqual(nucSeq[1].Count_Invalid_woIndel, 0);
            Assert.AreEqual(nucSeq[1].Count_A, 2);
            Assert.AreEqual(nucSeq[1].Count_T, 5);
            Assert.AreEqual(nucSeq[1].Count_C, 5);
            Assert.AreEqual(nucSeq[1].Count_G, 5);
            Assert.AreEqual(nucSeq[1].Count_N, 7);
            Assert.Contains("hallo", nucSeq[1].Names);
            Assert.Contains("hallo2", nucSeq[1].Names);
            Assert.AreEqual(nucSeq[1].Comment, "this is a comment");
            Assert.AreEqual(nucSeq[1].ParentName, "");
            Assert.AreEqual(nucSeq[1].SequenceCount, 1);




        }

        [Test]
        public void Test_RandomGenerator_NucleotideSequence()
        {
            ISequenceReader<NucleotideSequence> reader = new RandomGenerater_NucleotideSequence(10, 200, 20);
            List<NucleotideSequence> seq = reader.GetSequences(10);

            Assert.AreEqual(seq.Count, 10);
            Assert.AreEqual(seq[0].GC_Percent, 20.0);
            Assert.AreEqual(seq[9].GC_Percent, 20);
            Assert.AreEqual(seq[0].Length, 200);
            Assert.AreEqual(seq[9].Length, 200);

            reader = new RandomGenerater_NucleotideSequence(10, 150, 0);
            int count = 0;
            seq.Clear();
            List<NucleotideSequence> temp;
            for (int i = 0; i < 10; i++)
            {
                seq.Add(reader.GetNextSequence());

            }

  
            Assert.AreEqual(seq.Count, 10);
            Assert.AreEqual(seq[0].Length, 150);
            Assert.AreEqual(seq[4].Length, 150);
            Assert.AreEqual(seq[0].GC_Percent, 0);
            Assert.AreEqual(seq[9].GC_Percent, 0);

            reader = new RandomGenerater_NucleotideSequence(1, 150, 100);
            seq.Clear();
            seq = reader.GetSequences(1);

            Assert.AreEqual(seq.Count, 1);
            Assert.AreEqual(seq[0].Length, 150);
            Assert.AreEqual(seq[0].GC_Percent, 100);


            reader = new RandomGenerater_NucleotideSequence(1, 300, 80);
            seq.Clear();
            seq = reader.GetSequences(1);

            Assert.AreEqual(seq.Count, 1);
            Assert.AreEqual(seq[0].Length, 300);
            Assert.AreEqual(seq[0].GC_Percent, 80);



        }

        [Test]
        public void Test_SequenceReaderDecorator_ExtractFlankingSequence()
        {
            string path = @"..\..\TestFiles\Nunit_seq_NucleotideSequenceDecorator_ExtractFlankingRegion.txt";
            List<INucleotideSequence> seq = new List<INucleotideSequence>(4);
            seq.Add(new NucleotideSequence(null, "extractOne", new Sequence_ByteArray("CGCGCGCGCGCGCGCGCGCG"), 21, 40));
            seq.Add(new NucleotideSequence(null, "extractSecond", new Sequence_ByteArray("GATGATGATGATGAT"), 21, 35));
            seq.Add(new NucleotideSequence(null, "extractThird", new Sequence_ByteArray("TTTTTTTTTT"), 41, 50));
            seq.Add(new NucleotideSequence(null, "extractThird", new Sequence_ByteArray("CCCCCCCCCC"), 91, 100));

            //Bow in awe and look at the following code which converts the list into another list of the desired property
            //fucking difficult to understand, convertAll method needs a converter. The converter a delegate, I use a fucking ... delegate because i do not want to create a method here
            List<IPositionable> seqp = seq.ConvertAll(new Converter<INucleotideSequence, IPositionable>(delegate(INucleotideSequence nuc) { return (IPositionable)nuc; }));
            Dictionary<string, List<IPositionable>> pos = SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>.GroupPositionablesWithEqualParent(seqp);
            Assert.AreEqual(pos.Count, 3);
            Assert.IsTrue(pos.ContainsKey("extractOne"));
            Assert.IsTrue(pos.ContainsKey("extractSecond"));
            Assert.IsTrue(pos.ContainsKey("extractThird"));
            FastaReader_Sequence<NucleotideSequence> fr = new FastaReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCGN());
            ISequenceReader<NucleotideSequence> nr = new SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>(fr, seqp,false, 10, 5);

            NucleotideSequence temp;
            List<NucleotideSequence> s = new List<NucleotideSequence>();
            int count = 0;

            while ((temp = nr.GetNextSequence())!=null)
            {
                s.Add(temp);
                count++;
            }

            Assert.AreEqual(count, 4);
            Assert.AreEqual(s.Count, 4);

            Assert.AreEqual(s[0].Sequence.ToString(), "AAAAAAAAAACGCGCGCGCGCGCGCGCGCGAAAAA");
            Assert.AreEqual(s[0].Start, 11);
            Assert.AreEqual(s[0].End, 45);
            Assert.AreEqual(s[0].Length, 35);
            Assert.AreEqual(s[0].ParentName, "extractOne");
            Assert.AreEqual(s[0].Name, "extractOne");


            Assert.AreEqual(s[1].Sequence.ToString(), "TTTTTTTTTTGATGATGATGATGATTTTTT");
            Assert.AreEqual(s[1].Start, 11);
            Assert.AreEqual(s[1].End, 40);
            Assert.AreEqual(s[1].Length, 30);
            Assert.AreEqual(s[1].ParentName, "extractSecond");

            Assert.AreEqual(s[2].Sequence.ToString(), "AAAAAAAAAATTTTTTTTTTAAAAA");
            Assert.AreEqual(s[2].Start, 31);
            Assert.AreEqual(s[2].End, 55);
            Assert.AreEqual(s[2].Length, 25);
            Assert.AreEqual(s[2].ParentName, "extractThird");
            Assert.AreEqual(s[2].Name, "extractThird");

            Assert.AreEqual(s[3].Sequence.ToString(), "TTTTTTTTTTCCCCCCCCCCTTTTT");
            Assert.AreEqual(s[3].Start, 81);
            Assert.AreEqual(s[3].End, 105);
            Assert.AreEqual(s[3].Length, 25);
            Assert.AreEqual(s[3].ParentName, "extractThird");


            //Round TWO
            fr = new FastaReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCGN());
            nr = new SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>(fr, seqp, 40);

     
            s.Clear();
            s = nr.GetAllSequences();
            Assert.AreEqual(s.Count, 4);

            Assert.AreEqual(s[0].Sequence.ToString(), "AAAAAAAAAACGCGCGCGCGCGCGCGCGCGAAAAAAAAAA");
            Assert.AreEqual(s[0].Start, 11);
            Assert.AreEqual(s[0].End, 50);
            Assert.AreEqual(s[0].Length, 40);
            Assert.AreEqual(s[0].ParentName, "extractOne");

            Assert.AreEqual(s[1].Sequence.ToString(), "CCTTTTTTTTTTGATGATGATGATGATTTTTTTTTTTCCC");
            Assert.AreEqual(s[1].Start, 9);
            Assert.AreEqual(s[1].End, 48);
            Assert.AreEqual(s[1].Length, 40);
            Assert.AreEqual(s[1].ParentName, "extractSecond");

            Assert.AreEqual(s[2].Sequence.ToString(), "AAAAAAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAAAA");
            Assert.AreEqual(s[2].Start, 26);
            Assert.AreEqual(s[2].End, 65);
            Assert.AreEqual(s[2].Length, 40);
            Assert.AreEqual(s[2].ParentName, "extractThird");

            Assert.AreEqual(s[3].Sequence.ToString(), "TTTTTTTTTTTTTTTCCCCCCCCCCTTTTTTTTTTTTTTT");
            Assert.AreEqual(s[3].Start, 76);
            Assert.AreEqual(s[3].End, 115);
            Assert.AreEqual(s[3].Length, 40);
            Assert.AreEqual(s[3].ParentName, "extractThird");



            //Round Three 
            fr = new FastaReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCGN());
            nr = new SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>(fr, seqp, 30, 0);

            s.Clear();

            s = nr.GetAllSequences();

            Assert.AreEqual(s.Count, 4);

            Assert.AreEqual(s[0].Sequence.ToString(), "TTTTTTTTTTAAAAAAAAAACGCGCGCGCGCGCGCGCGCGAAAAAAAAAA");
            Assert.AreEqual(s[0].Start, 1);
            Assert.AreEqual(s[0].End, 50);
            Assert.AreEqual(s[0].Length, 50);
            Assert.AreEqual(s[0].ParentName, "extractOne");

            //Round Four 
            fr = new FastaReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCGN());
            nr = new SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>(fr, seqp, 400);


            List<NucleotideSequence> tempList;
            s.Clear();
            count = 0;
            while ((tempList = nr.GetSequences(2)).Count > 0)
            {
                s.AddRange(tempList);
                tempList.Clear();
                count++;

            }
            Assert.AreEqual(count, 2);
            Assert.AreEqual(s.Count, 4);

            Assert.AreEqual(s[0].Sequence.ToString(), "TTTTTTTTTTAAAAAAAAAACGCGCGCGCGCGCGCGCGCGAAAAAAAAAATTTTTTTTTT");
            Assert.AreEqual(s[0].Start, 1);
            Assert.AreEqual(s[0].End, 60);
            Assert.AreEqual(s[0].Length, 60);
            Assert.AreEqual(s[0].ParentName, "extractOne");

            Assert.AreEqual(s[3].Sequence.ToString(), "TTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAA");
            Assert.AreEqual(s[3].Start, 1);
            Assert.AreEqual(s[3].End, 140);
            Assert.AreEqual(s[3].Length, 140);
            Assert.AreEqual(s[3].ParentName, "extractThird");


            //Round Five 
            fr = new FastaReader_Sequence<NucleotideSequence>(path, new DNASequenceValidator_ATCGN());
            nr = new SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>(fr, seqp, true, 400);


            s.Clear();

            s = nr.GetAllSequences();
            Assert.AreEqual(s.Count, 4);

            Assert.AreEqual(s[0].Sequence.ToString(), "TTTTTTTTTTAAAAAAAAAACGCGCGCGCGCGCGCGCGCGAAAAAAAAAATTTTTTTTTT");
            Assert.AreEqual(s[0].Start, 1);
            Assert.AreEqual(s[0].End, 60);
            Assert.AreEqual(s[0].Length, 60);
            Assert.AreEqual(s[0].ParentName, "extractOne");


            Assert.AreEqual(s[3].Sequence.ToString(), "TTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAATTTTTTTTTTAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTCCCCCCCCCCTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAA");
            Assert.AreEqual(s[3].Start, 1);
            Assert.AreEqual(s[3].End, 140);
            Assert.AreEqual(s[3].Length, 140);
            Assert.AreEqual(s[3].ParentName, "extractThird");

        }

        [Test]
        public void Test_SequenceReaderDecorater_UniqueSequenceID()
        {
            string path = @"..\..\TestFiles\Nunit_seq_NucleotideSequenceDecorator_UniqueSeqID.txt";
            ISequenceReader<NucleotideSequence> sr = new FastaReader_Sequence<NucleotideSequence>(path, new SequenceValidator_AllValid());
            sr = new SequenceReaderDecorator_UniqueSeqID<NucleotideSequence>(sr);
            NucleotideSequence temp;
            List<NucleotideSequence> seq = new List<NucleotideSequence>();
            int count = 0;
            while ((temp = sr.GetNextSequence())!=null)
            {
                seq.Add(temp);
                count++;
            }
            Assert.AreEqual(count, 10);
            Assert.AreEqual(seq.Count, 10);
            CollectionAssert.AllItemsAreUnique(seq);
            Dictionary<string, bool> contains = new Dictionary<string, bool>(seq.Count * 2);
            for (int i = 0; i < seq.Count; i++)
            {
                Assert.IsFalse(contains.ContainsKey(seq[i].Name));
                contains.Add(seq[i].Name, true);
                Console.WriteLine("Name number {0}: {1}", i + 1, seq[i].Name);
            }

            //Round two
            sr = new FastaReader_Sequence<NucleotideSequence>(path, new SequenceValidator_AllValid());
            sr = new SequenceReaderDecorator_UniqueSeqID<NucleotideSequence>(sr);
            seq.Clear();
            seq = sr.GetAllSequences();
     
            Assert.AreEqual(seq.Count, 10);
            CollectionAssert.AllItemsAreUnique(seq);
            contains = new Dictionary<string, bool>(seq.Count * 2);
            for (int i = 0; i < seq.Count; i++)
            {
                Assert.IsFalse(contains.ContainsKey(seq[i].Name));
                contains.Add(seq[i].Name, true);
                Console.WriteLine("Name number {0}: {1}", i + 1, seq[i].Name);
            }
        }

        [Test]
        public void Test_FastasReader_NucleotideSequence()
        {

            string[] paths = new string[5];
            paths[0] = @"..\..\TestFiles\Nunit_seq_fastreader_nucleotideSequence.txt";
            paths[1] = @"..\..\TestFiles\Nunit_seq_fastreader_nucseqInfo.txt";
            paths[2] = @"..\..\TestFiles\Nunit_seq_fastaReaderChopper_nucleotideSequence.txt";
            paths[3] = @"..\..\TestFiles\Nunit_seq_fastawriter_nucleotideSequence.txt";
            paths[4] = @"..\..\TestFiles\Nunit_seqMisc_fastasreader_nucleotideSequence.txt";
            FastaMultiFileReader_Sequence<NucleotideSequence> frNuc = new FastaMultiFileReader_Sequence<NucleotideSequence>(paths, new SequenceValidator_AllValid());
            frNuc.SequencePrototype = new Sequence_String();

            List<NucleotideSequence> seq = new List<NucleotideSequence>();
            NucleotideSequence temp;
            int count = 0;
            while ((temp = frNuc.GetNextSequence())!=null)
            {
                count++;
                seq.Add(temp);
            }
            frNuc.Close();

            Assert.AreEqual(count, 19);
            Assert.AreEqual(seq.Count, 19);
            Assert.AreEqual(seq[15].Name, "sixteen");
            Assert.AreEqual(seq[16].Name, "seventeen");
            Assert.AreEqual(seq[17].Name, "eighteen");
            Assert.AreEqual(seq[18].Name, "nineteen");
            Assert.IsTrue(seq[0].Sequence is Sequence_String);
            Assert.IsTrue(seq[18].Sequence is Sequence_String);


            ///Round two
            seq.Clear();
            frNuc = null;
            frNuc = new FastaMultiFileReader_Sequence<NucleotideSequence>(paths, new SequenceValidator_AllValid());
            frNuc.SequencePrototype = new Sequence_ByteArray();
            seq = frNuc.GetAllSequences();
            Assert.AreEqual(seq.Count, 19);
            Assert.AreEqual(seq[15].Name, "sixteen");
            Assert.AreEqual(seq[16].Name, "seventeen");
            Assert.AreEqual(seq[17].Name, "eighteen");
            Assert.AreEqual(seq[18].Name, "nineteen");
            Assert.IsTrue(seq[0].Sequence is Sequence_ByteArray);
            Assert.IsTrue(seq[18].Sequence is Sequence_ByteArray);


            //Round Three
            paths = new string[1];
            paths[0] = @"G:\Programs\00TestSequencen\Nunit_seq_fastreader_nucleotideSequence.txt";
            frNuc = new FastaMultiFileReader_Sequence<NucleotideSequence>(paths, new DNASequenceValidator_ATCGNYRWSKMBDHV(),FastaReaderOptions.IgnoreInvalid);


            seq.Clear();
            seq = frNuc.GetAllSequences();
            Assert.AreEqual(seq.Count, 2);
            Assert.AreEqual(seq[0].Name, "first");
            Assert.AreEqual(seq[0].Sequence.ToString(), "AAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGNNNNNNNNNNNNNNNNNNNNNAAAAAAAAAATTTTTTTTTTGGGGGGGGGGCCCCCCCCCC");

            Assert.AreEqual(seq[1].Name, "third");
            Assert.AreEqual(seq[1].Sequence.ToString(), "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");




        }
       
    }

    [TestFixture]
    public class Seq_Misc
    {

        public static void Start()
        {
            Seq_Misc p = new Seq_Misc();
            p.Test_CreateUniqueSeqID();
            p.Test_Trimming();
        }

        [Test]
        public void Test_CreateUniqueSeqID()
        {
            string[] names = new string[] { "franz", "franz", "franz", "hans", "hans", "hans", "samson", "helga", "martin" };
            List<string> namesList = new List<string>(names.Length);


            foreach (string name in names)
            {
                namesList.Add(name);
            }

            Dictionary<string,System.Collections.Generic.Queue<string>> matrix = CreateUniqueSeqID.GetUniqueIDMatrix(namesList);
            Assert.AreEqual(matrix.Keys.Count, 5);

            CreateUniqueSeqID seqId = new CreateUniqueSeqID(matrix);

            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_1");
            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_2");
            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_3");
            Assert.AreEqual(seqId.GetUniqueID("samson"), "samson");
            Assert.AreEqual(seqId.GetUniqueID("martin"), "martin");
            Assert.AreEqual(seqId.GetUniqueID("helga"), "helga");

            namesList.Add("franz_1");
            namesList.Add("franz_2");
            namesList.Add("franz_3");
            namesList.Add("franz_1*");

            matrix = CreateUniqueSeqID.GetUniqueIDMatrix(namesList);
            Assert.AreEqual(matrix.Keys.Count, 9);
            seqId = new CreateUniqueSeqID(matrix);
            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_1**");
            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_2*");
            Assert.AreEqual(seqId.GetUniqueID("franz"), "franz_3*");
            Assert.AreEqual(seqId.GetUniqueID("franz_1"), "franz_1");
            Assert.AreEqual(seqId.GetUniqueID("franz_1*"), "franz_1*");


        }

        [Test]
        public void Test_Trimming()
        {
            NucleotideSequence ns;
            NucleotideSequence trimed;




            ns = new NucleotideSequence("fritz", new Sequence_ByteArray("ACGTTTTTAAAAAGCA"));
            ns.QualitySequence = new QualitySequence("fritz", new Sequence_ByteArray("\a\b\t\n\v\f\r"));

            trimed = TrimNucleotideSequence.Trim_5Prime_and_3Prime(ns, 3, 4);
            Assert.AreEqual(trimed.Name, "fritz");
            Assert.AreEqual(trimed.QualitySequence.Name, "fritz");
            Assert.AreEqual(trimed.Sequence.ToString(), "TTTTTAAAA");
            Assert.AreEqual(trimed.QualitySequence.Sequence.ToString(), "\a\b\t\n\v\f");

            trimed = TrimNucleotideSequence.Trim_5Prime_and_3Prime(ns, 0, 0);
            Assert.AreEqual(trimed.Name, "fritz");
            Assert.AreEqual(trimed.QualitySequence.Name, "fritz");
            Assert.AreEqual(trimed.Sequence.ToString(), "ACGTTTTTAAAAAGCA");
            Assert.AreEqual(trimed.QualitySequence.Sequence.ToString(), "\a\b\t\n\v\f\r");

            ns = new NucleotideSequence("fritz", new Sequence_ByteArray("ACGTTTTTTTTTTTTTTTAACAAAAAGCAT"));//30 length
            ns.QualitySequence = new QualitySequence("fritz", new Sequence_ByteArray("\a\b\t\n\v\f\r"));

            trimed = TrimNucleotideSequence.Trim_5Prime_and_3PrimePolyA(ns, 3, 50, 5);
            Assert.AreEqual(trimed.Name, "fritz");
            Assert.AreEqual(trimed.QualitySequence.Name, "fritz");
            Assert.AreEqual(trimed.Sequence.ToString(), "TTTTTTTTTTTTTTT");
            Assert.AreEqual(trimed.QualitySequence.Sequence.ToString(), "\a\b\t\n\v\f\r");


            trimed = TrimNucleotideSequence.Trim_5Prime_and_3PrimePolyA(ns, 0, 50, 6);
            Assert.AreEqual(trimed.Name, "fritz");
            Assert.AreEqual(trimed.QualitySequence.Name, "fritz");
            Assert.AreEqual(trimed.Sequence.ToString(), "ACGTTTTTTTTTTTTTTTAACAAAAAGCAT");
            Assert.AreEqual(trimed.QualitySequence.Sequence.ToString(), "\a\b\t\n\v\f\r");




        }


    }

    [TestFixture]
    public class Seq_Sort
    {
        public static void Start()
        {
            Seq_Sort p = new Seq_Sort();
            p.Test_Sort_ParentName_Start_Ascending();
            p.Test_Sort_ParentName_Start_Descending();
            p.Test_Sort_Length_Ascending();
            p.Test_Sort_Length_Descending();
        }

        [Test]
        public void Test_Sort_ParentName_Start_Ascending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 300, 312));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 55, 60));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 1, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 4, 3));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 30, 312));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 15, 60));

            nucS.Sort(new SortPositionables_ParentName_StartPos_Ascending<NucleotideSequence>());
            Assert.AreEqual(nucS[0].ParentName, "franz");
            Assert.AreEqual(nucS[0].Start, 12);
            Assert.AreEqual(nucS[1].ParentName, "franz");
            Assert.AreEqual(nucS[1].Start, 55);
            Assert.AreEqual(nucS[2].ParentName, "franz");
            Assert.AreEqual(nucS[2].Start, 111);
            Assert.AreEqual(nucS[3].ParentName, "franz");
            Assert.AreEqual(nucS[3].Start, 300);

            
            Assert.AreEqual(nucS[4].ParentName, "frenz");
            Assert.AreEqual(nucS[4].Start, 12);
            Assert.AreEqual(nucS[5].ParentName, "frenz");
            Assert.AreEqual(nucS[5].Start, 15);
            Assert.AreEqual(nucS[6].ParentName, "frenz");
            Assert.AreEqual(nucS[6].Start, 30);
            Assert.AreEqual(nucS[7].ParentName, "frenz");
            Assert.AreEqual(nucS[7].Start, 111);

            Assert.AreEqual(nucS[8].ParentName, "xavier");
            Assert.AreEqual(nucS[8].Start, 1);
            Assert.AreEqual(nucS[9].ParentName, "xavier");
            Assert.AreEqual(nucS[9].Start, 2);
            Assert.AreEqual(nucS[10].ParentName, "xavier");
            Assert.AreEqual(nucS[10].Start, 3);
            Assert.AreEqual(nucS[11].ParentName, "xavier");
            Assert.AreEqual(nucS[11].Start, 4);
        }

        [Test]
        public void Test_Sort_ParentName_Start_Descending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 300, 312));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 55, 60));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 1, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 4, 3));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 30, 312));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 15, 60));


            nucS.Sort(new SortPositionables_ParentName_StartPos_Descending<NucleotideSequence>());
            Assert.AreEqual(nucS[11].ParentName, "franz");
            Assert.AreEqual(nucS[11].Start, 12);
            Assert.AreEqual(nucS[10].ParentName, "franz");
            Assert.AreEqual(nucS[10].Start, 55);
            Assert.AreEqual(nucS[9].ParentName, "franz");
            Assert.AreEqual(nucS[9].Start, 111);
            Assert.AreEqual(nucS[8].ParentName, "franz");
            Assert.AreEqual(nucS[8].Start, 300);


            Assert.AreEqual(nucS[7].ParentName, "frenz");
            Assert.AreEqual(nucS[7].Start, 12);
            Assert.AreEqual(nucS[6].ParentName, "frenz");
            Assert.AreEqual(nucS[6].Start, 15);
            Assert.AreEqual(nucS[5].ParentName, "frenz");
            Assert.AreEqual(nucS[5].Start, 30);
            Assert.AreEqual(nucS[4].ParentName, "frenz");
            Assert.AreEqual(nucS[4].Start, 111);

            Assert.AreEqual(nucS[3].ParentName, "xavier");
            Assert.AreEqual(nucS[3].Start, 1);
            Assert.AreEqual(nucS[2].ParentName, "xavier");
            Assert.AreEqual(nucS[2].Start, 2);
            Assert.AreEqual(nucS[1].ParentName, "xavier");
            Assert.AreEqual(nucS[1].Start, 3);
            Assert.AreEqual(nucS[0].ParentName, "xavier");
            Assert.AreEqual(nucS[0].Start, 4);
        }

        [Test]
        public void Test_Sort_Length_Ascending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACCG"), 1, 7));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACAACC"), 1, 9));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATAATAC"), 1, 8));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACC"), 1, 6));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 1, 4));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("AT"), 1, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATCC"), 1, 5));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("A"), 1, 2));

            nucS.Sort(new SortPositionables_Length_Ascending<NucleotideSequence>());
            Assert.AreEqual(nucS[0].ParentName,"xavier");
            Assert.AreEqual(nucS[0].End, 2);
            Assert.AreEqual(nucS[1].ParentName, "xavier");
            Assert.AreEqual(nucS[1].End, 3);
            Assert.AreEqual(nucS[2].ParentName, "xavier");
            Assert.AreEqual(nucS[2].End, 4);
            Assert.AreEqual(nucS[3].ParentName, "xavier");
            Assert.AreEqual(nucS[3].End, 5);

            Assert.AreEqual(nucS[4].ParentName, "franz");
            Assert.AreEqual(nucS[4].End, 6);
            Assert.AreEqual(nucS[5].ParentName, "franz");
            Assert.AreEqual(nucS[5].End, 7);
            Assert.AreEqual(nucS[6].ParentName, "franz");
            Assert.AreEqual(nucS[6].End, 8);
            Assert.AreEqual(nucS[7].ParentName, "franz");
            Assert.AreEqual(nucS[7].End, 9);
        }

        [Test]
        public void Test_Sort_Length_Descending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACCG"), 1, 7));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACAACC"), 1, 9));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATAATAC"), 1, 8));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATACC"), 1, 6));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 1, 4));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("AT"), 1, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATCC"), 1, 5));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("A"), 1, 2));

            nucS.Sort(new SortPositionables_Length_Descending<NucleotideSequence>());
            Assert.AreEqual(nucS[7].ParentName, "xavier");
            Assert.AreEqual(nucS[7].End, 2);
            Assert.AreEqual(nucS[6].ParentName, "xavier");
            Assert.AreEqual(nucS[6].End, 3);
            Assert.AreEqual(nucS[5].ParentName, "xavier");
            Assert.AreEqual(nucS[5].End, 4);
            Assert.AreEqual(nucS[4].ParentName, "xavier");
            Assert.AreEqual(nucS[4].End, 5);

            Assert.AreEqual(nucS[3].ParentName, "franz");
            Assert.AreEqual(nucS[3].End, 6);
            Assert.AreEqual(nucS[2].ParentName, "franz");
            Assert.AreEqual(nucS[2].End, 7);
            Assert.AreEqual(nucS[1].ParentName, "franz");
            Assert.AreEqual(nucS[1].End, 8);
            Assert.AreEqual(nucS[0].ParentName, "franz");
            Assert.AreEqual(nucS[0].End, 9);
        }


        /*

        [Test]
        public void Test_Sort_ParentName_Start_IndelShift_Ascending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 15, 2, 0));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 312, 3, 1));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 111, 115, 1, 1));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 55, 60, 1, 1));

            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3, 2, 2));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3, 3, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3, 1, 1));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3, 3, 3));

            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 30, 312));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 15, 60));

            nucS.Sort(new SortPositionables_ParentName_StartPos_IndelShift_Ascending<NucleotideSequence>());
            Assert.AreEqual(nucS[0].ParentName, "franz");
            Assert.AreEqual(nucS[0].Start, 12);
            Assert.AreEqual(nucS[0].Start_IndelShift, 2);
            Assert.AreEqual(nucS[1].ParentName, "franz");
            Assert.AreEqual(nucS[1].Start, 12);
            Assert.AreEqual(nucS[1].Start_IndelShift, 3);
            Assert.AreEqual(nucS[2].ParentName, "franz");
            Assert.AreEqual(nucS[2].Start, 55);
            Assert.AreEqual(nucS[2].Start_IndelShift, 1);
            Assert.AreEqual(nucS[3].ParentName, "franz");
            Assert.AreEqual(nucS[3].Start, 111);
            Assert.AreEqual(nucS[3].Start_IndelShift, 1);


            Assert.AreEqual(nucS[4].ParentName, "frenz");
            Assert.AreEqual(nucS[4].Start, 12);
            Assert.AreEqual(nucS[4].Start_IndelShift, null);
            Assert.AreEqual(nucS[5].ParentName, "frenz");
            Assert.AreEqual(nucS[5].Start, 15);
            Assert.AreEqual(nucS[6].ParentName, "frenz");
            Assert.AreEqual(nucS[6].Start, 30);
            Assert.AreEqual(nucS[7].ParentName, "frenz");
            Assert.AreEqual(nucS[7].Start, 111);
            Assert.AreEqual(nucS[7].Start_IndelShift, null);

            Assert.AreEqual(nucS[8].ParentName, "xavier");
            Assert.AreEqual(nucS[8].Start, 2);
            Assert.AreEqual(nucS[8].Start_IndelShift, 2);
            Assert.AreEqual(nucS[9].ParentName, "xavier");
            Assert.AreEqual(nucS[9].Start, 2);
            Assert.AreEqual(nucS[9].Start_IndelShift, 3);
            Assert.AreEqual(nucS[10].ParentName, "xavier");
            Assert.AreEqual(nucS[10].Start, 3);
            Assert.AreEqual(nucS[10].Start_IndelShift, 1);
            Assert.AreEqual(nucS[11].ParentName, "xavier");
            Assert.AreEqual(nucS[11].Start, 3);
            Assert.AreEqual(nucS[11].Start_IndelShift, 3);
        }

        [Test]
        public void Test_Sort_ParentName_Start_IndelShift_Descending()
        {
            List<NucleotideSequence> nucS = new List<NucleotideSequence>();
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 15, 2, 0));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 12, 312, 3, 1));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 111, 115, 1, 1));
            nucS.Add(new NucleotideSequence("", "franz", new Sequence_ByteArray("ATA"), 55, 60, 1, 1));

            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3, 2, 2));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 2, 3, 3, 3));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3, 1, 1));
            nucS.Add(new NucleotideSequence("", "xavier", new Sequence_ByteArray("ATA"), 3, 3, 3, 3));

            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 12, 15));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 30, 312));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 111, 115));
            nucS.Add(new NucleotideSequence("", "frenz", new Sequence_ByteArray("ATA"), 15, 60));

            nucS.Sort(new SortPositionables_ParentName_StartPos_IndelShift_Descending<NucleotideSequence>());
            Assert.AreEqual(nucS[11].ParentName, "franz");
            Assert.AreEqual(nucS[11].Start, 12);
            Assert.AreEqual(nucS[11].Start_IndelShift, 2);
            Assert.AreEqual(nucS[10].ParentName, "franz");
            Assert.AreEqual(nucS[10].Start, 12);
            Assert.AreEqual(nucS[10].Start_IndelShift, 3);
            Assert.AreEqual(nucS[9].ParentName, "franz");
            Assert.AreEqual(nucS[9].Start, 55);
            Assert.AreEqual(nucS[9].Start_IndelShift, 1);
            Assert.AreEqual(nucS[8].ParentName, "franz");
            Assert.AreEqual(nucS[8].Start, 111);
            Assert.AreEqual(nucS[8].Start_IndelShift, 1);


            Assert.AreEqual(nucS[7].ParentName, "frenz");
            Assert.AreEqual(nucS[7].Start, 12);
            Assert.AreEqual(nucS[7].Start_IndelShift, null);
            Assert.AreEqual(nucS[6].ParentName, "frenz");
            Assert.AreEqual(nucS[6].Start, 15);
            Assert.AreEqual(nucS[5].ParentName, "frenz");
            Assert.AreEqual(nucS[5].Start, 30);
            Assert.AreEqual(nucS[4].ParentName, "frenz");
            Assert.AreEqual(nucS[4].Start, 111);
            Assert.AreEqual(nucS[4].Start_IndelShift, null);

            Assert.AreEqual(nucS[3].ParentName, "xavier");
            Assert.AreEqual(nucS[3].Start, 2);
            Assert.AreEqual(nucS[3].Start_IndelShift, 2);
            Assert.AreEqual(nucS[2].ParentName, "xavier");
            Assert.AreEqual(nucS[2].Start, 2);
            Assert.AreEqual(nucS[2].Start_IndelShift, 3);
            Assert.AreEqual(nucS[1].ParentName, "xavier");
            Assert.AreEqual(nucS[1].Start, 3);
            Assert.AreEqual(nucS[1].Start_IndelShift, 1);
            Assert.AreEqual(nucS[0].ParentName, "xavier");
            Assert.AreEqual(nucS[0].Start, 3);
            Assert.AreEqual(nucS[0].Start_IndelShift, 3);
        }

        */




    }


    [TestFixture]
    public class Seq_Statistics
    {
        public static void Start()
        {
            Seq_Statistics p = new Seq_Statistics();
            p.Test_ChromosomeBin();
            p.Test_ChromosomeDistributionStatisitc();
            p.Test_MannWhitneyU();
            p.Test_PoissonDistribution();
            p.Test_SequenceStatistic();
        }


        [Test]
        public void Test_SequenceStatistic()
        {
            List<NucleotideSequence> nucs = new List<NucleotideSequence>();
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTTT"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTT"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("CCCCCCCCCCC"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("GGGGG"), 12, 1));

            NucleotideSequenceStatistic<NucleotideSequence> stat = new NucleotideSequenceStatistic<NucleotideSequence>(nucs);
            Assert.AreEqual(stat.Count, 4);
            Assert.AreEqual(stat.AverageGC_Content, 0.5);
            Assert.AreEqual(stat.AverageLength, 8);
            Assert.AreEqual(stat.TotalLengthOfAllNucleotides, 32);
            Assert.AreEqual(stat.StdDev_Length_RandomVariable, 3.0);
            Assert.AreEqual(stat.NucSeqInfo, null);
            Assert.AreEqual(stat.DensityAllNucleotides, null);
            Assert.AreEqual(stat.DensityATUCG, null);
            Assert.AreEqual(stat.DensityWithoutN, null);
            Assert.AreEqual(stat.FilesPerCount, null);
            Assert.AreEqual(stat.NucleotideSequences, nucs);
            Assert.AreEqual(stat.NucleotidesPerCount, null);
            Assert.AreEqual(stat.StdDev_Length_SampleOfPopulation, Math.Sqrt(12.0));

            NucSeqInfo ns = new NucSeqInfo("franz", 100000000, 10000000);
            ns.Length = 1000;
            ns.Count_C = 100;
            ns.Count_G = 100;
            ns.Count_A = 100;
            ns.Count_U = 100;
            ns.Count_T = 100;
            ns.Count_N = 200;
            ns.Count_Special = 300;

            stat = new NucleotideSequenceStatistic<NucleotideSequence>(nucs, ns);

            Assert.AreEqual(stat.Count, 4);
            Assert.AreEqual(stat.AverageGC_Content, 0.5);
            Assert.AreEqual(stat.AverageLength, 8);
            Assert.AreEqual(stat.TotalLengthOfAllNucleotides, 32);
            Assert.AreEqual(stat.StdDev_Length_RandomVariable, 3.0);
            Assert.AreEqual(stat.StdDev_Length_SampleOfPopulation, Math.Sqrt(12.0));

            Assert.AreEqual(stat.NucSeqInfo, ns);
            Assert.AreEqual(stat.DensityAllNucleotides, 4000);
            Assert.AreEqual(stat.DensityATUCG, 8000);
            Assert.AreEqual(stat.DensityWithoutN, 5000);
            Assert.AreEqual(stat.FilesPerCount, 0.25);
            Assert.AreEqual(stat.NucleotideSequences, nucs);
            Assert.AreEqual(stat.NucleotidesPerCount, 250);






        }

        [Test]
        public void Test_ChromosomeBin()
        {
            List<NucleotideSequence> nucs = new List<NucleotideSequence>();
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTTT"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTT"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("CCCCCCCCCCC"), 12, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("GGGGG"), 12, 1));

            ChromosomeBin<NucleotideSequence> stat = new ChromosomeBin<NucleotideSequence>(nucs, "hans", 1000, 1000);
            Assert.AreEqual(stat.Count, 4);

            Assert.AreEqual(stat.AverageLength, 8);
            Assert.AreEqual(stat.TotalLengthOfAllNucleotides, 32);
            Assert.AreEqual(stat.StdDev_Length_RandomVariable, 3.0);

            Assert.AreEqual(stat.DensityAllNucleotides, 4000);
            Assert.AreEqual(stat.DensityATUCG, null);
            Assert.AreEqual(stat.DensityWithoutN, null);
            Assert.AreEqual(stat.FilesPerCount, 0.25);
            Assert.AreEqual(stat.NucleotideSequences, nucs);
            Assert.AreEqual(stat.NucleotidesPerCount, 250);
            Assert.AreEqual(stat.StdDev_Length_SampleOfPopulation, Math.Sqrt(12.0));
            Assert.AreEqual(stat.Start, 1000);
            Assert.AreEqual(stat.End, 1999);
            Assert.AreEqual(stat.Length, 1000);
            Assert.AreEqual(stat.PositionMiddle, 1499.5);

        }


        [Test]
        public void Test_ChromosomeDistributionStatisitc()
        {
            List<NucleotideSequence> nucs = new List<NucleotideSequence>();
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTTT"), 1, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTT"), 250, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("CCCCCCCCCCC"), 750, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("GGGGG"), 1000, 1));

            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTTTTTTTT"), 1001, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTT"), 1250, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("CCCCCCCCCCCCCCCC"), 1750, 1));
            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("GGGGGGGGGG"), 2000, 1));

            nucs.Add(new NucleotideSequence("leng10", "franz", new Sequence_ByteArray("TTTTTTTTTTTGGGGGGGGGGG"), 2001, 2010));
            nucs.Add(new NucleotideSequence("leng10", "ferdinand", new Sequence_ByteArray("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"), 0, 1));
            nucs.Add(new NucleotideSequence("leng10", "fritz", new Sequence_ByteArray("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"), 500, 1));
            nucs.Add(new NucleotideSequence("leng10", "fritz", new Sequence_ByteArray("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"), 1500, 1));

            ChromosomeStatisticProvider<NucleotideSequence> provider = new ChromosomeStatisticProvider<NucleotideSequence>(nucs);
            ChromosomeDistributionStatistic<NucleotideSequence> stat = provider.GetBinStatistic("franz", 1000, 1000);
            Assert.AreEqual(stat.NucleotideSequences.Count, 9);
            Assert.AreEqual(stat.Count, 3);
            Assert.AreEqual(stat[0].AverageLength, 8);
            Assert.AreEqual(stat[0].Count, 4);
            Assert.AreEqual(stat[0].Start, 1);
            Assert.AreEqual(stat[0].End, 1000);

            Assert.AreEqual(stat[1].Count, 4);
            Assert.AreEqual(stat[1].Start, 1001);
            Assert.AreEqual(stat[1].End, 2000);
            Assert.AreEqual(stat[1].AverageLength, 13);
            Assert.AreEqual(stat[1].StdDev_Length_RandomVariable, 3.0);


            Assert.AreEqual(stat[2].Count, 1);
            Assert.AreEqual(stat[2].AverageLength, 22);

            //
            //ROUND TWO
            //
            stat = provider.GetBinStatistic("franz", 1000, 300);

            Assert.AreEqual(stat.NucleotideSequences.Count, 9);
            Assert.AreEqual(stat.Count, 5);

            Assert.AreEqual(stat[0].AverageLength, 8);
            Assert.AreEqual(stat[0].Count, 4);

            Assert.AreEqual(stat[0].Start, 1);
            Assert.AreEqual(stat[0].End, 1000);

            Assert.AreEqual(stat[1].Count, 4);
            Assert.AreEqual(stat[2].Count, 4);
            Assert.AreEqual(stat[2].Start, 601);
            Assert.AreEqual(stat[2].End, 1600);
            Assert.AreEqual(stat[3].Count, 4);


            Assert.AreEqual(stat[4].End, 2010);
            Assert.AreEqual(stat[4].Start, 1201);
            Assert.AreEqual(stat[4].Count, 4);

            //
            //ROUND THREE
            //

            stat = provider.GetBinStatistic("fritz", 2000, 1000);

            Assert.AreEqual(stat.Count, 1);
            Assert.AreEqual(stat.NucleotideSequences.Count, 2);
        }



        [Test]
        public void Test_MannWhitneyU()
        {
            List<double> basis = new List<double>();
            basis.Add(10); basis.Add(12); basis.Add(14);
            List<double> left = new List<double>();
            left.Add(9); left.Add(11); left.Add(13);
            List<double> right = new List<double>();
            right.Add(11); right.Add(13); right.Add(15);

            MannWhitneyU bl = new MannWhitneyU(basis, left);
            MannWhitneyU br = new MannWhitneyU(basis, right);
            MannWhitneyU bb = new MannWhitneyU(basis, basis);
            Assert.AreEqual(bl.PositionOfToTest, -1);
            Assert.AreEqual(br.PositionOfToTest, 1);

            Assert.AreEqual(bb.PositionOfToTest, 0);
            Assert.AreEqual(bb.P_OneSided, 0.5);
            Assert.AreEqual(bb.P_TwoSided, 1.0);

            ///
            ///STEP2
            ///
            List<double> kontrollgruppe = new List<double>();
            List<double> versuchsgruppe = new List<double>();
            kontrollgruppe.Add(74.0); kontrollgruppe.Add(68.0); kontrollgruppe.Add(53); kontrollgruppe.Add(94); kontrollgruppe.Add(80);
            kontrollgruppe.Add(54);
            Assert.AreEqual(kontrollgruppe.Count, 6);
            versuchsgruppe.Add(86); versuchsgruppe.Add(55); versuchsgruppe.Add(77); versuchsgruppe.Add(92);
            versuchsgruppe.Add(63); versuchsgruppe.Add(89); versuchsgruppe.Add(90); versuchsgruppe.Add(93);
            Assert.AreEqual(versuchsgruppe.Count, 8);
            MannWhitneyU buch = new MannWhitneyU(kontrollgruppe, versuchsgruppe);

            ///This values are obtained assuming a large sample size, this is called normal approximation
            ///Values are calculated from http://elegans.swmed.edu/~leon/stats/utest.html
            Assert.Less(Math.Abs(buch.P_OneSided - 0.122639), 0.000001);
            Assert.Less(Math.Abs(buch.P_TwoSided - 0.245278), 0.000001);
            Assert.AreEqual(buch.PositionOfToTest, 1);
        }

        [Test]
        public void Test_PoissonDistribution()
        {
            PoissonDistribution distri = new PoissonDistribution(12.5);

            Assert.IsTrue(Math.Abs(distri.CumulativeProbability_Right(13) - 0.372164656) < 0.0000001);
            Assert.IsTrue(Math.Abs(distri.CumulativeProbability_Right(14) - 0.274968116) < 0.000001);

            Assert.IsTrue(Math.Abs(distri.CumulativeProbability_Left(13) - 0.627835344) < 0.00000001);
            Assert.IsTrue(Math.Abs(distri.CumulativeProbability_Left(14) - 0.725031884) < 0.00000001);

            Assert.IsTrue(Math.Abs(distri.Probability(13) - 0.108860125) < 0.00000001);
            Assert.IsTrue(Math.Abs(distri.Probability(14) - 0.09719654) < 0.00000001);

        }


    }
}