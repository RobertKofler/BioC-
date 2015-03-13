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
using System.IO;

namespace Bio.Seq.IO
{

    [Flags]
    public enum FastaReaderOptions
    {
        None=0x00,
        TrimExcess = 0x01,
        IgnoreInvalid = 0x02,
        ConvertToUppercase=0x04
    }

    [Flags]
    public enum FastaReaderNucSeqInfoOptions
    {
        None=0x00,
        IgnoreInvalid = 0x01,
        FullDetails = 0x02
    }


#if full
    /// <summary>
    /// Reading fasta files, return a collection of INucleotideSequences;
    /// Comments are also read (;)
    /// A ISequenceValidatior may be passed
    /// </summary>
    public class FastaReader_Sequence<TSimpleSequence> : ISequenceReader<TSimpleSequence>
        where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
    {
        private StreamReader sr;
        private readonly ISequenceValidator validator;
        private readonly bool ignoreInvalid = false;
        private readonly bool trimExcess = false;
        private string name;
        private TSimpleSequence templateSimpleSequence = new TSimpleSequence();

        /// <summary>
        /// seqTemplate is used as template for sequences. Acts like a factory, the NewSequence() mehtod is called to obtain a new ISequence of the same type
        /// Therefore the user does not have to know the type of the Sequence
        /// </summary>
        private ISequenceContainer seqTemplate = SequenceFactory.GetDefaultSequence();

        private List<TSimpleSequence> sequences;


        public FastaReader_Sequence(string path, ISequenceValidator validator, bool ignoreInvalid, int bufferSize) : this(new StreamReader(path, Encoding.ASCII), validator, ignoreInvalid, bufferSize) { }
        public FastaReader_Sequence(StreamReader sr, ISequenceValidator validator, bool ignoreInvalid, int bufferSize)
        {
            this.sr = sr;
            this.validator = validator;
            this.bufferSize = bufferSize;
            this.ignoreInvalid = ignoreInvalid;
            name = "";
            sequences = new List<TSimpleSequence>(bufferSize);
        }

        private void ParseOptions(FastaReaderOptions options)
        {
            if (options & FastaReaderOptions.IgnoreInvalid == FastaReaderOptions.IgnoreInvalid) this.ignoreInvalid = true;
            if (options & FastaReaderOptions.TrimExcess == FastaReaderOptions.TrimExcess) this.trimExcess = true;
        }




        public FastaReader_Sequence(string path, ISequenceValidator validator, bool ignoreInvalid)
            : this(new StreamReader(path, Encoding.ASCII), validator, ignoreInvalid)
        {
        }
        public FastaReader_Sequence(StreamReader sr, ISequenceValidator validator, bool ignoreInvalid)
        {
            name = "";
            this.sr = sr;
            this.validator = validator;
            this.ignoreInvalid = ignoreInvalid;
            sequences = new List<TSimpleSequence>(this.bufferSize);
        }

        /// <summary>
        /// Returns a number of sequences as specified through the buffer size.
        /// If the size of a sequence is greater than >200000 than the buffer size is set to 10
        /// if the size of the sequence is greater than 2000000 the buffer size for this file is set to 1
        /// </summary>
        /// <returns></returns>
        public List<TSimpleSequence> GetSequences()
        {
#if DEBUG
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
            sw.Start();
#endif
            //The Mark updated
            sequences.Clear();
            string line = "";
            string comment = null;
            ISequenceContainer sequence = seqTemplate.NewSequence();

            //As long as there are lines in rhe file proceed!
            while ((sequences.Count < bufferSize) &&
                ((line = sr.ReadLine()) != null))
            {
                //Get rid of space
                line = line.Trim();

                //If line is empty or starts with the comment symbol skip the line
                if (line == "") continue;

                //The > resets
                if (line.StartsWith(">"))
                {
                    //Store old FastaFile
                    if (name != "" && sequence.Length != 0)
                    {
                        //Adjust bufferSize for large sequences
                        if (sequence.Length > 2000000 && bufferSize > 1) bufferSize = 1;
                        else if (sequence.Length > 200000 && bufferSize > 10) bufferSize = 10;

                        if (validator.IsValid(sequence))
                        {
                            sequence.TrimExcess();
                            sequences.Add(templateSimpleSequence.NewInstance(name, comment, sequence));
                        }
                        else if (!ignoreInvalid) throw new InvalidSequenceException(validator);
                    }

                    //Create new name and sequence
                    name = line.Substring(1, line.Length - 1);
                    sequence = seqTemplate.NewSequence();
                    comment = null;
                }
                else if (line.StartsWith(";"))
                {
                    comment = comment ?? ""; //Comment is the content of comment except comment=null then comment=""
                    comment += line.Substring(1, line.Length - 1);
                }
                else
                {
                    sequence.Append(line);
                }

            }
            //Also store last file
            if (name != "" && sequence.Length != 0)
            {
                if (validator.IsValid(sequence))
                {
                    sequence.TrimExcess();
                    sequences.Add(templateSimpleSequence.NewInstance(name, comment, sequence));
                }
                else if (!ignoreInvalid) throw new InvalidSequenceException(validator);
            }

#if DEBUG
            sw.Stop();
            System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReader_NucleotideSequence; {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, seqTemplate.GetType().ToString()));


#endif


            return sequences;
        }

        public TSimpleSequence GetNextSequence()
        {
        }
        public List<TSimpleSequence> GetAllSequences()
        {
        }

        public ISequenceContainer SequencePrototype
        {
            get
            {
                return this.seqTemplate;
            }
            set
            {
                this.seqTemplate = value;
            }
        }

        public void Close()
        {
            this.sr.Close();
        }

    }

#endif

    /// <summary>
    /// Interface for a class reading ISimpleSequences; Generic class
    /// SimpleSequence readers may either read one sequence at the time or all at once
    /// </summary>
    /// <typeparam name="TSimpleSequence"></typeparam>
    public interface ISequenceWriter<TSimpleSequence> where TSimpleSequence: ISimpleSequence
    {
          void WriteSequences(List<TSimpleSequence> sequences);
          void WriteSequence(TSimpleSequence sequence);
          void Close();
    }


     /// <summary>
     /// The following interface has to be implemented by SimpleSequence derived classes, if they should be read with a ISequenceReader
     /// Implementation of this interface allows them to be read with ISequenceReader 
     /// </summary>
     /// <typeparam name="TSimpleSequence"></typeparam>
     public interface IFastaFileReadable<TSimpleSequence> where TSimpleSequence : ISimpleSequence
     {
         TSimpleSequence NewInstance(string name, string comment, ISequenceContainer sequence);
         TSimpleSequence NewInstance(string name, string parentName, ISequenceContainer sequence, int startPosition, int endPosition);


     }

     public interface INucSeqInfoReader
     {
         List<INucSeqInfo> GetSequences(int number);
         List<INucSeqInfo> GetAllSequences();
         INucSeqInfo GetNextSequence();
         void Close();
     }
     public interface INucSeqInfoWriter
     {
         void WriteSeqInfos(List<INucSeqInfo> sequenceInfos);
         void WriteSeqInfo(INucSeqInfo sequenceInfo);
         void Close();
     }


     public interface ISequenceReader<TSimpleSequence> where TSimpleSequence: ISimpleSequence
     {
         /// <summary>
         /// Retrieve all sequences at once
         /// </summary>
         /// <returns></returns>
         List<TSimpleSequence> GetAllSequences();

         /// <summary>
         /// Retrieve only the next sequence
         /// </summary>
         /// <returns></returns>
         TSimpleSequence GetNextSequence();

         /// <summary>
         /// Retrieve the specified number of sequences
         /// </summary>
         /// <param name="number"></param>
         /// <returns></returns>
         List<TSimpleSequence> GetSequences(int number);


         ISequenceContainer SequencePrototype { get;set;}
         void Close();

     
     }

     /// <summary>
     /// Reading fasta files, return a collection of INucleotideSequences;
     /// Comments are also read (;)
     /// A ISequenceValidatior may be passed
     /// </summary>
     public class FastaReader_Sequence<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
     {
         private StreamReader sr;
         private readonly ISequenceValidator validator;
         private readonly bool ignoreInvalid = false;
         private readonly bool trimExcess = false;
         private readonly bool convertToUppercase = false;
         private string name;
         private TSimpleSequence templateSimpleSequence = new TSimpleSequence();

         /// <summary>
         /// seqTemplate is used as template for sequences. Acts like a factory, the NewSequence() mehtod is called to obtain a new ISequence of the same type
         /// Therefore the user does not have to know the type of the Sequence
         /// </summary>
         private ISequenceContainer seqTemplate = SequenceFactory.GetDefaultSequence();


         public FastaReader_Sequence(string path, ISequenceValidator validator,FastaReaderOptions options) : this(new StreamReader(path, Encoding.ASCII), validator, options) { }


         public FastaReader_Sequence(StreamReader sr, ISequenceValidator validator, FastaReaderOptions options)
         {
             this.sr = sr;
             this.validator = validator;

             if ((options & FastaReaderOptions.IgnoreInvalid) == FastaReaderOptions.IgnoreInvalid) ignoreInvalid = true;
             if ((options & FastaReaderOptions.TrimExcess) == FastaReaderOptions.TrimExcess) trimExcess = true;
             if ((options & FastaReaderOptions.ConvertToUppercase)==FastaReaderOptions.ConvertToUppercase )convertToUppercase=true;

             name = "";
         }




         public FastaReader_Sequence(string path, ISequenceValidator validator)
             : this(new StreamReader(path, Encoding.ASCII), validator)
         {
         }
         public FastaReader_Sequence(StreamReader sr, ISequenceValidator validator)
         {
             name = "";
             this.sr = sr;
             this.validator = validator;
         }

         /// <summary>
         /// Returns the next sequence
         /// </summary>
         /// <returns></returns>
         public TSimpleSequence GetNextSequence()
         {

             string line = "";
             string comment = null;
             bool invalid = true;
             ISequenceContainer sequence = seqTemplate.NewSequence();

             //As long as there are lines in rhe file proceed!
             while ( (line = sr.ReadLine()) != null)
             {
                 //Get rid of space
                 line = line.Trim();

                 //If line is empty or starts with the comment symbol skip the line
                 if (line == "") continue;

                 //The > resets
                 if (line.StartsWith(">"))
                 {
                     //Store old FastaFile
                     if (name != "" && sequence.Length != 0)
                     {

                         if (validator.IsValid(sequence))
                         {
                             if(trimExcess) sequence.TrimExcess();
                             if (convertToUppercase) sequence.ToUpper();
                             return templateSimpleSequence.NewInstance(name, comment, sequence);

                         }
                         else if (!ignoreInvalid)
                         {
                             throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);
                         }
                     }

                     //Create new name and sequence
                     name = line.Substring(1, line.Length - 1);
                     sequence = seqTemplate.NewSequence();
                     comment = null;
                 }
                 else if (line.StartsWith("#"))
                 {
                     comment = comment ?? ""; //Comment is the content of comment except comment=null then comment=""
                     comment += line.Substring(1, line.Length - 1);
                 }
                 else
                 {
                     sequence.Append(line);
                 }

             }

             //Also store last file
             if (name != "" && sequence.Length != 0)
             {
                 if (validator.IsValid(sequence))
                 {
                     if(trimExcess) sequence.TrimExcess();
                     if (convertToUppercase) sequence.ToUpper();
                     TSimpleSequence seq= templateSimpleSequence.NewInstance(name, comment, sequence);
                     name = "";
                     sequence = null;
                     return seq;
                 }
                 else if (!ignoreInvalid)
                 {
                     throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);
                 }
             }


             return default(TSimpleSequence);
         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return this.seqTemplate;
             }
             set
             {
                 this.seqTemplate = value;
             }
         }

         public void Close()
         {
             this.sr.Close();
         }



         public List<TSimpleSequence> GetAllSequences()
         {
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             List<TSimpleSequence> ret = new List<TSimpleSequence>();
             TSimpleSequence next;
             while ((next = GetNextSequence()) != null)
             {
                 ret.Add(next);
             }
             ret.TrimExcess();

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReader; All {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, ret.Count, sw.ElapsedMilliseconds, ret.Count <= 0 ? 0 : ret[0].Length, seqTemplate.GetType().ToString()));
#endif
             return ret;
         }



         public List<TSimpleSequence> GetSequences(int number)
         {

             List<TSimpleSequence> ret = new List<TSimpleSequence>();
             TSimpleSequence next;
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             while (ret.Count<number && (next = GetNextSequence()) != null)
             {
                 ret.Add(next);
             }
             ret.TrimExcess();

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReader; A batch of {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, ret.Count, sw.ElapsedMilliseconds, ret.Count <= 0 ? 0 : ret[0].Length, seqTemplate.GetType().ToString()));
#endif
             return ret;
         }


     }


     /// <summary>
     /// Writes a simple nucleotide Sequence in a fasta file, heading the recommendation that the one line should have a length of 
     /// at most 50 nucleotides
     /// </summary>
     public class FastaWriter<TSimpleSequence> : ISequenceWriter<TSimpleSequence> where TSimpleSequence : ISimpleSequence
     {
         private StreamWriter sw;
         public FastaWriter(string path) : this(new StreamWriter(path, false, Encoding.ASCII)) { }
         public FastaWriter(StreamWriter sw)
         {
             this.sw = sw;
         }

         /// <summary>
         /// Writes nucleotide sequences into a fasta-style output file.
         /// Since the fasta definition also contains a comment line (starts with ;) the comments are also writen to the output file
         /// If the .Comment property is null or "" no comment line will be written
         /// </summary>
         /// <param name="sequences"></param>
         public void WriteSequences(List<TSimpleSequence> sequences)
         {
#if DEBUG
             System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
             stopw.Start();
#endif
             foreach (TSimpleSequence ns in sequences)
             {
                 WriteSequence(ns);
             }

#if DEBUG
             stopw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaWriter; {1} files saved in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, stopw.ElapsedMilliseconds,
                 sequences.Count <= 0 ? 0 : sequences[0].Length, sequences.Count <= 0 ? "no sequence" : sequences[0].Sequence.GetType().ToString()));
#endif
         }

         public void Close()
         {
             this.sw.Close();
         }



      


         public void WriteSequence(TSimpleSequence sequence)
         {
             int count = 0;
             // string line = "";
             ISequenceContainer seq = sequence.Sequence;
             sw.WriteLine('>' + sequence.Name);
             if (sequence.Comment != null && sequence.Comment != "") sw.WriteLine('#' + sequence.Comment);

             while (seq.Length >= count + 50)
             {
                 sw.WriteLine(seq.ToString(count, 50));
                 count += 50;
             }
             sw.WriteLine(seq.ToString(count, seq.Length - count));
             sw.WriteLine();
         }


     }

     /// <summary>
     /// Writes a simple nucleotide Sequence into a byte-file;
     /// Attention the byte(255) and byte(254) are used for internal puposes and mark the end of the sequence and the end of the sequence-name respectively
     /// the first 4000 bytes are reserved for the name of the sequence
     /// the second 4000 bytes are reserved for comments
     /// </summary>
     public class ByteWriter<TSimpleSequence> : ISequenceWriter<TSimpleSequence> where TSimpleSequence : ISimpleSequence
     {
         private FileStream fs;
         public ByteWriter(string path) : this(new FileStream(path, FileMode.Create, FileAccess.Write)) { }
         public ByteWriter(FileStream fs)
         {
             this.fs = fs;
         }

         /// <summary>
         /// Writes nucleotide sequences into a byte-output file.
         /// byte(255) is reserved and marks the end of a sequence 
         /// first 4000 bytes are for the name the next 4000 bytes for the comment
         /// </summary>
         /// <param name="sequences"></param>
         public void WriteSequences(List<TSimpleSequence> sequences)
         {
#if DEBUG
             System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
             stopw.Start();
#endif
             foreach (TSimpleSequence ns in sequences)
             {
                 WriteSequence(ns);

             }

#if DEBUG
             stopw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} ByteWriter; {1} files saved in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, stopw.ElapsedMilliseconds,
                 sequences.Count <= 0 ? 0 : sequences[0].Length, sequences.Count <= 0 ? "no sequence" : sequences[0].Sequence.GetType().ToString()));


#endif
         }


         public void WriteSequence(TSimpleSequence sequence)
         {

             string seQName = sequence.Name;
             string seqComment = sequence.Comment;

             byte[] name = new byte[seQName.Length + 1];
             byte[] comment = new byte[1];
             if (seqComment != null) comment = new byte[seqComment.Length + 1];

             //First write name of the sequence
             for (int i = 0; i < seQName.Length; i++)
             {
#if DEBUG
                 if ((((byte)seQName[i]) == (byte)255) || (((byte)seQName[i]) == (byte)254)) throw new Exception("crap");
#endif
                 name[i] = (byte)seQName[i];

             }
             name[name.Length - 1] = (byte)254;

             //Second write the comment
             for (int i = 0; i < (comment.Length - 1); i++)
             {
#if DEBUG
                 if ((((byte)seqComment[i]) == (byte)255) || (((byte)seqComment[i]) == (byte)254)) throw new Exception("crap");
#endif
                 comment[i] = (byte)seqComment[i];
             }
             comment[comment.Length - 1] = (byte)254;

             fs.Write(name, 0, name.Length);
             fs.Write(comment, 0, comment.Length);

             //Third write the sequence to the file
             ISequenceContainer seq = sequence.Sequence;
             for (int i = 0; i < seq.Length; i++)
             {
#if DEBUG
                 if ((((byte)seq[i]) == (byte)255) || (((byte)seq[i]) == (byte)254)) throw new Exception("crap");
#endif

                 fs.WriteByte((byte)seq[i]);
             }
             fs.Write(new byte[] { (byte)255 }, 0, 1);
         }

         public void Close()
         {
             this.fs.Close();
         }
     }


     /// <summary>
     /// Reading byte files, return a collection of INucleotideSequences;
     /// byte(255) denotes the end of a sequence
     /// A ISequenceValidatior may be passed
     /// </summary>
     public class ByteReader_Sequence<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
     {
         private FileStream fs;
         private readonly ISequenceValidator validator;
         private readonly bool ignoreInvalid = false;
         private readonly bool convertToUppercase = false;
         private readonly bool trimExcess = false;

         private TSimpleSequence templateSimpleSequence = new TSimpleSequence();


         /// <summary>
         /// seqTemplate is used as template for sequences. Acts like a factory, the NewSequence() mehtod is called to obtain a new ISequence of the same type
         /// Therefore the user does not have to know the type of the Sequence
         /// </summary>
         private ISequenceContainer seqTemplate = SequenceFactory.GetDefaultSequence();

         private List<TSimpleSequence> sequences;


         public ByteReader_Sequence(string path, ISequenceValidator validator, FastaReaderOptions options) : this(new FileStream(path, FileMode.Open, FileAccess.Read), validator, options) { }


         public ByteReader_Sequence(FileStream fs, ISequenceValidator validator, FastaReaderOptions options)
         {
             this.fs = fs;
             this.validator = validator;

             if ((options & FastaReaderOptions.IgnoreInvalid) != 0) ignoreInvalid = true;
             if ((options & FastaReaderOptions.TrimExcess) != 0) trimExcess = true;
             if ((options & FastaReaderOptions.ConvertToUppercase) != 0) convertToUppercase = true;


         }




         public ByteReader_Sequence(string path, ISequenceValidator validator)
             : this(new FileStream(path, FileMode.Open, FileAccess.Read), validator)
         {
         }


         public ByteReader_Sequence(FileStream fs, ISequenceValidator validator)
         {
             this.fs = fs;
             this.validator = validator;
         }



         public TSimpleSequence GetNextSequence()
         {


             byte by;
             string name = null;
             string comment = null;
             bool invalid = true;
             ISequenceContainer sequence = seqTemplate.NewSequence();
             bool readName = true;
             TSimpleSequence ret = default(TSimpleSequence);



             //As long as the buffer can still hold sequences
             while (invalid)
             {
                 if (readName) //If this is a new sequence -> read the first 8000 bytes
                 {

                     StringBuilder sb = new StringBuilder("");

                     //The end of the file is a int with the value -1
                     //can only be obtained if one byte is read
                     int peek = fs.ReadByte();
                     int byt = 0;
                     if (peek == -1)
                     {
                         break;

                     }
                     else if (peek == (byte)254)
                     {
                         name = "";
                     }
                     else //if the obtained character was not the end of the file insert it into the byte array
                     {
                         sb.Append((char)peek);
                         while ((byt = fs.ReadByte()) != (byte)254)
                         {
                             sb.Append((char)byt);
                         }


                     }
                     name = sb.ToString();
                     sb = new StringBuilder("");

                     while ((byt = fs.ReadByte()) != (byte)254)
                     {
                         sb.Append((char)byt);
                     }
                     if (sb.Length > 0) comment = sb.ToString();

                     readName = false;


                 }
                 else
                 {
                     by = (byte)fs.ReadByte();

                     if (by == (byte)255) //if the read character ist the sequence terminating character store the sequence 
                     {

                         if (validator.IsValid(sequence))
                         {
                             if (trimExcess) sequence.TrimExcess();
                             if (convertToUppercase) sequence.ToUpper();

                             ret=templateSimpleSequence.NewInstance(name, comment, sequence);
                             invalid = false;
                         }
                         else if (!ignoreInvalid) throw new InvalidSequenceException(validator, validator.GetInvalid(sequence),name);


                         name = null;
                         comment = null;
                         readName = true;
                         sequence = seqTemplate.NewSequence();

                     }
                     else sequence.Append(by); //if the read byte is not sequence terminating append it to the sequence

                 }
             }
             return ret;
         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return this.seqTemplate;
             }
             set
             {
                 this.seqTemplate = value;
             }
         }





         public void Close()
         {
             this.fs.Close();
         }



  

         public List<TSimpleSequence> GetAllSequences()
         {
             List<TSimpleSequence> ret = new List<TSimpleSequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             TSimpleSequence seq;
             while ((seq = GetNextSequence()) != null)
             {
                 ret.Add(seq);
             }
             ret.TrimExcess();


#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} ByteReader_Sequence; All {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, ret.Count, sw.ElapsedMilliseconds, ret.Count <= 0 ? 0 : ret[0].Length, seqTemplate.GetType().ToString()));
#endif

             return ret;

         }


         public List<TSimpleSequence> GetSequences(int number)
         {
             List<TSimpleSequence> ret = new List<TSimpleSequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             TSimpleSequence seq;
             while (ret.Count < number && (seq = GetNextSequence()) != null)
             {
                 ret.Add(seq);
             }
             ret.TrimExcess();


#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} ByteReader_Sequence; A batch of {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, ret.Count, sw.ElapsedMilliseconds, ret.Count <= 0 ? 0 : ret[0].Length, seqTemplate.GetType().ToString()));
#endif

             return ret;
         }

  
     }

    


     /// <summary>
     /// Reading fasta files, return a collection of INucSeqInfo;
     /// This class only reads general informations for nucleotide sequences contained in a fasta file like: are the characters uppercase; how many A,T,C,G etc
     /// Comments are also read (;) 
     /// The detail level may be specified; fullDetails = false -> only the fastaID and the comments are returned
     /// A ISequenceValidatior may be passed
     /// </summary>
     public class FastaReader_NucSeqInfo : INucSeqInfoReader
     {
         private StreamReader sr;
         private readonly ISequenceValidator validator;


         INucSeqInfo nsi = new NucSeqInfo("");
         private bool fullDetails = false;
         private readonly bool ignoreInvalid = false;


         /// <summary>
         /// Constructor for reading only the fasta ID, the length and the fasta comment,
         /// all fasta sequences are valid
         /// </summary>
         public FastaReader_NucSeqInfo(StreamReader sr, ISequenceValidator validator)
         {
             this.sr = sr;
             this.validator = validator;
         }

         /// <summary>
         /// Constructor for reading only the fasta ID, the length and the fasta comment,
         /// all fasta sequences are valid
         /// </summary>
         public FastaReader_NucSeqInfo(string path,ISequenceValidator validator)
         {
             sr = new StreamReader(path, Encoding.ASCII);
             this.validator = validator;
         }

         /// <summary>
         /// Constructor for reading only the fasta ID, the length and the fasta comment,
         /// all fasta sequences are valid
         /// </summary>
         public FastaReader_NucSeqInfo(string path)
             : this(path, new SequenceValidator_AllValid())
         {
 
         }

         /// <summary>
         /// Constructor for reading only the fasta ID, the length and the fasta comment,
         /// all fasta sequences are valid
         /// </summary>
         public FastaReader_NucSeqInfo(StreamReader sr)
             : this(sr, new SequenceValidator_AllValid())
         {
         }





         public FastaReader_NucSeqInfo(string path, ISequenceValidator validator, FastaReaderNucSeqInfoOptions options)
             : this(new StreamReader(path, Encoding.ASCII), validator, options) { }


         public FastaReader_NucSeqInfo(StreamReader sr, ISequenceValidator validator, FastaReaderNucSeqInfoOptions options)
             : this()
         {
             this.sr = sr;
             this.validator = validator;

             if ((options & FastaReaderNucSeqInfoOptions.FullDetails) ==FastaReaderNucSeqInfoOptions.FullDetails) this.fullDetails = true;
             if ((options & FastaReaderNucSeqInfoOptions.IgnoreInvalid) == FastaReaderNucSeqInfoOptions.IgnoreInvalid) this.ignoreInvalid = true;
         }


         private FastaReader_NucSeqInfo()
         {

         }






         public INucSeqInfo GetNextSequence()
         {

             string line = "";
             // string name = "";
             string comment = null;
             bool invalid = true;
             INucSeqInfo ret = null;


             //As long as there are lines in rhe file proceed!
             while (invalid && (line = sr.ReadLine()) != null)
             {
                 //Get rid of space
                 line = line.Trim();

                 //If line is empty or starts with the comment symbol skip the line
                 if (line == "") continue;

                 //The > resets
                 if (line.StartsWith(">"))
                 {
                     //Store old FastaFile
                     if (nsi.Name != "" && nsi.Length != 0) //
                     {
                         nsi.Comment = comment;
                         if (validator.IsValid(nsi))
                         {
                             ret=nsi.Copy();
                             invalid = false;
                         }
                         else if (!ignoreInvalid) throw new InvalidSequenceException(validator,new List<int>(),nsi.Name);


                     }
                     nsi = new NucSeqInfo(line.Substring(1, line.Length - 1));
                     nsi.IsUpper = true;
                     comment = "";

                 }
                 else if (line.StartsWith("#"))
                 {
                     comment = comment ?? ""; //Comment is the content of comment except comment=null then comment=""
                     comment += line.Substring(1, line.Length - 1);
                 }
                 else
                 {
                     //In any case increase count the length
                     nsi.Length += line.Length;


                     if (fullDetails)
                     {
                         for (int i = 0; i < line.Length; i++)
                         {
                             char c = line[i];
                             if (Char.IsLower(c)) nsi.IsUpper = false;

                             switch (c)
                             {


                                 case 'a':
                                 case 'A': nsi.Count_A++;
                                     break;
                                 case 't':
                                 case 'T': nsi.Count_T++;
                                     break;
                                 case 'c':
                                 case 'C': nsi.Count_C++;
                                     break;
                                 case 'g':
                                 case 'G': nsi.Count_G++;
                                     break;
                                 case 'u':
                                 case 'U': nsi.Count_U++;
                                     break;
                                 case 'n':
                                 case 'N': nsi.Count_N++;
                                     break;
                                 case '-': nsi.Count_Indel++;
                                     break;



                                 //Special symbols:

                                 case 'Y':   //Pyrmidin (C or T)
                                 case 'y':
                                 case 'r':   //Purin (G or A)
                                 case 'R':
                                 case 'w':   //Weich (A or T)
                                 case 'W':
                                 case 'S':   //Stark(C or G)
                                 case 's':
                                 case 'k':   //Ketogruppe (G or T)
                                 case 'K':
                                 case 'm':   //Aminogruppe (A or C)
                                 case 'M':
                                 case 'B':   //Nicht A (C,G,T)
                                 case 'b':
                                 case 'd':   //Nicht C (A,T,G)
                                 case 'D':
                                 case 'h':   //Nicht G (A,T,C)
                                 case 'H':
                                 case 'v':   //Nicht T (A,C,G)
                                 case 'V': nsi.Count_Special++;
                                     break;
                                 default: nsi.Count_Invalid_woIndel++;
                                     break;



                             }

                         }
                     }
                     else
                     {
                         nsi.IsUpper = null;
                     }

                 }

             }
             //Also store last file
             if (nsi.Name != "" && nsi.Length != 0)//
             {
                 nsi.Comment = comment;
                 if (validator.IsValid(nsi))
                 {
                    ret=nsi.Copy();
                 }
                 else if (!ignoreInvalid) throw new InvalidSequenceException(validator,new List<int>(),nsi.Name);
                 nsi = new NucSeqInfo("");
             }




             return ret;
         }

         public void Close()
         {
             this.sr.Close();
         }






         public List<INucSeqInfo> GetSequences(int number)
         {
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             List<INucSeqInfo> sequences = new List<INucSeqInfo>();
             INucSeqInfo next;
             while (sequences.Count < number && (next = this.GetNextSequence()) != null)
             {
                 sequences.Add(next);
             }

             

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReader_NucSeqInfo; A batch of {1} files loaded in {2} ms; First file contained {3} nucleotides; Full detail level: {4}", DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, fullDetails));

#endif
             return sequences;
         }



         public List<INucSeqInfo> GetAllSequences()
         {
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif

             List<INucSeqInfo> sequences = new List<INucSeqInfo>();
             INucSeqInfo next;

             while ((next = this.GetNextSequence()) != null)
             {
                 sequences.Add(next);
             }


#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReader_NucSeqInfo; All {1} files loaded in {2} ms; First file contained {3} nucleotides; Full detail level: {4}", DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, fullDetails));

#endif

             return sequences;
         }




     }


    public class FastaMultiFileReader_Sequence<TSimpleSequence> : ISequenceReader<TSimpleSequence> where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
    {
        private ISequenceContainer template = SequenceFactory.GetDefaultSequence();
        private ISequenceValidator validator;
        private FastaReaderOptions options;
        private StreamReader[] srs;
        private int fileNumber = 0;
        private FastaReader_Sequence<TSimpleSequence> activeReader;


        public FastaMultiFileReader_Sequence(string[] paths, ISequenceValidator validator, FastaReaderOptions options)
        {
            if (paths.Length < 1) throw new InvalidDataException("paths (string[] does not contain a valid path");
            srs = new StreamReader[paths.Length];
            for (int i = 0; i < paths.Length; i++)
            {
                srs[i] = new StreamReader(paths[i], Encoding.ASCII);
            }
            this.validator = validator;
            this.options = options;

            this.NewReader();
        }

        public FastaMultiFileReader_Sequence(string[] paths, ISequenceValidator validator)
        {
            if (paths.Length < 1) throw new InvalidDataException("paths (string[] does not contain a valid path");

            srs = new StreamReader[paths.Length];
            for (int i = 0; i < paths.Length; i++)
            {
                srs[i] = new StreamReader(paths[i], Encoding.ASCII);
            }
            this.validator = validator;


            this.NewReader();
        }


        public TSimpleSequence GetNextSequence()
        {

            TSimpleSequence temp;
            while ((temp = activeReader.GetNextSequence())== null)
            {
                srs[fileNumber].Close(); //Close old StreamReader
                fileNumber++;

                if (fileNumber < srs.Length)
                {
                    NewReader();
                }
                else
                {
                    return default(TSimpleSequence);
                }
            }
            return temp;
        }



        private void NewReader()
        {
            activeReader = new FastaReader_Sequence<TSimpleSequence>(srs[fileNumber], validator, options);
            this.activeReader.SequencePrototype = this.template;
        }

        public ISequenceContainer SequencePrototype
        {
            get
            {
                return this.template;
            }
            set
            {
                this.template = value;
                activeReader.SequencePrototype = value;
            }
        }

        public void Close()
        {
            this.activeReader.Close();
        }

        public List<TSimpleSequence> GetAllSequences()
        {
            List<TSimpleSequence> toRet = new List<TSimpleSequence>();
            TSimpleSequence temp;
            while ((temp = this.GetNextSequence()) != null)
            {
                toRet.Add(temp);
            }
            return toRet;
        }


        public List<TSimpleSequence> GetSequences(int number)
        {
            List<TSimpleSequence> toRet = new List<TSimpleSequence>();
            TSimpleSequence temp;
            while (toRet.Count<number && (temp = this.GetNextSequence()) != null)
            {
                toRet.Add(temp);
            }
            return toRet;
        }
    }


     /// <summary>
     /// Reads nucleotide sequences, and chops the sequence into chunks of the given size.
     /// A overhead might be specified (i.e. a small nucleotide sequence that overlaps with the previous sequence).
     /// All chunks derived from a large sequence have the same fasta-ID and comment.
     /// </summary>
     public class FastaReaderSequenceChopper<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
     {
         private readonly ISequenceValidator validator;
         private readonly bool ignoreInvalid = false;
         private readonly bool trimExcess = false;
         private readonly bool convertToUppercase = false;

         private StreamReader sr;
         private ISequenceContainer seqTemplate;

         //Variables necessary when operating
         private string name;
         private string comment;
         private ISequenceContainer sequence;
         private string previousOverhead = "";

         //Special FastaReader_Chopper variables
         private readonly int sequenceSize;
         private readonly int overhead;
         private TSimpleSequence templateSimpleSequence = new TSimpleSequence();


         public FastaReaderSequenceChopper(string path, ISequenceValidator validator, FastaReaderOptions options, int sequenceSize, int overhead)
             :
             this(new StreamReader(path, Encoding.ASCII), validator, options, sequenceSize, overhead) { }

         public FastaReaderSequenceChopper(StreamReader sr, ISequenceValidator validator, FastaReaderOptions options, int sequenceSize, int overhead)
             : this()
         {
             this.sr = sr;
             this.validator = validator;

             //setting operating varibles to their default values
             name = "";

             //Setting special variables
             if (overhead >= sequenceSize) throw new ArgumentOutOfRangeException("Overhead > Sequence size");
             this.sequenceSize = sequenceSize;
             this.overhead = overhead;

             if ((options & FastaReaderOptions.IgnoreInvalid) == FastaReaderOptions.IgnoreInvalid) ignoreInvalid = true;
             if ((options & FastaReaderOptions.TrimExcess) == FastaReaderOptions.TrimExcess) trimExcess = true;
             if ((options & FastaReaderOptions.ConvertToUppercase) == FastaReaderOptions.ConvertToUppercase) convertToUppercase = true;
         }



         private FastaReaderSequenceChopper()
         {
             this.seqTemplate = SequenceFactory.GetDefaultSequence();
             this.sequence = seqTemplate.NewSequence();
         }

         public FastaReaderSequenceChopper(string path, ISequenceValidator validator, FastaReaderOptions options)
             :
             this(new StreamReader(path, Encoding.ASCII), validator, options) { }


         public FastaReaderSequenceChopper(StreamReader sr, ISequenceValidator validator, FastaReaderOptions options)
             : this()
         {
             this.sr = sr;
             this.validator = validator;
             //setting operating varibles to their default values
             name = "";

             //
             this.sequenceSize = 50000000;
             this.overhead = 50;

             if ((options & FastaReaderOptions.IgnoreInvalid) == FastaReaderOptions.IgnoreInvalid) ignoreInvalid = true;
             if ((options & FastaReaderOptions.TrimExcess) == FastaReaderOptions.TrimExcess) trimExcess = true;
             if ((options & FastaReaderOptions.ConvertToUppercase) == FastaReaderOptions.ConvertToUppercase) convertToUppercase = true;

         }


         public TSimpleSequence GetNextSequence()
         {

             TSimpleSequence ret = default(TSimpleSequence);
             bool invalid = true;
             string line = "";

             //Append previous overhead
             sequence.Append(previousOverhead.Trim());


             //As long as there are lines in rhe file proceed!
             while (invalid &&  ((line = sr.ReadLine()) != null))
             {
                 //Get rid of space
                 line = line.Trim();
                 //If line is empty or starts with the comment symbol skip it
                 if (line == "") continue;

                 //The > resets
                 if (line.StartsWith(">"))
                 {
                     //Store old FastaFile
                     if (name != "" && !sequence.Equals(seqTemplate.NewSequence()))
                     {

                         if (validator.IsValid(sequence))
                         {
                             if (trimExcess) sequence.TrimExcess();
                             if (convertToUppercase) sequence.ToUpper();
                             ret= templateSimpleSequence.NewInstance(name, comment, sequence);
                             invalid = false;
                         }
                         else if (!ignoreInvalid) throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);
                     }

                     //Create new name
                     sequence = seqTemplate.NewSequence();
                     previousOverhead = "";
                     name = line.Substring(1, line.Length - 1);
                     comment = null;

                 }
                 else if (line.StartsWith("#"))
                 {
                     comment = comment ?? ""; //Comment is the content of comment except comment=null then comment=""
                     comment += line.Substring(1, line.Length - 1);
                 }
                 else if (sequence.Length >= this.sequenceSize)//If large files
                 {
                     //at first append the actual line
                     sequence.Append(line);

                     //get overhead for the next file eg sequence.Length=200052; sequenceSize=200000; overhead=50
                     //200052+50-200000=102; start position = 199999-50=199949
                     //10-5=5 
                     //01234      0123456789
                     //TTTTTAAAAA AAAAATTTTT 20+5-10=15
                     previousOverhead = sequence.ToString(sequenceSize - overhead, sequence.Length + overhead - sequenceSize);


                     if (validator.IsValid(sequence))
                     {
                         if (trimExcess) sequence.TrimExcess();
                         if (convertToUppercase) sequence.ToUpper();
                         ret=templateSimpleSequence.NewInstance(name, comment, sequence.SubSequence(0, sequenceSize));
                         invalid = false;
                     }
                     else if (!ignoreInvalid) throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);

                     sequence = seqTemplate.NewSequence();


                 }
                 else sequence.Append(line);

             }
             //Also store last file
             if (name != "" && !sequence.Equals(seqTemplate.NewSequence()))//ReadLine null
             {

                 if (validator.IsValid(sequence))
                 {
                     if(trimExcess) sequence.TrimExcess();
                     if(convertToUppercase) sequence.ToUpper();
                     ret = templateSimpleSequence.NewInstance(name, comment, sequence);
                 }
                 else if (!ignoreInvalid) throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);

                 previousOverhead = "";
                 sequence = seqTemplate.NewSequence();
             }



             return ret;
         }



         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return this.seqTemplate;
             }
             set
             {
                 this.seqTemplate = value;
                 this.sequence = seqTemplate.NewSequence();
             }
         }

         public void Close()
         {
             this.sr.Close();
         }




         public List<TSimpleSequence> GetAllSequences()
         {
             List<TSimpleSequence> sequences = new List<TSimpleSequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             TSimpleSequence temp;
             while ((temp = this.GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }



#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReaderChopper_NucleotideSequence; All {1} files loaded in {2} ms; First file contained {3} nucleotides; Specified sequnce size: {4}; Specified overhead: {5}; Sequence used: {6}",
                 DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, this.sequenceSize, overhead, seqTemplate.GetType().ToString()));

#endif
             return sequences;
        
         }



         public List<TSimpleSequence> GetSequences(int number)
         {
             List<TSimpleSequence> sequences = new List<TSimpleSequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             TSimpleSequence temp;
             while (sequences.Count<number && (temp = this.GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }



#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} FastaReaderChopper_NucleotideSequence; A batch of {1} files loaded in {2} ms; First file contained {3} nucleotides; Specified sequnce size: {4}; Specified overhead: {5}; Sequence used: {6}",
                 DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, this.sequenceSize, overhead, seqTemplate.GetType().ToString()));

#endif
             return sequences;
         }

 

     }

    [Flags]
    public enum QualityFileReaderOptions
    {
        None=0x00,
        IgnoreInvalid=0x01,
        TrimExcess=0x02
    }

    /// <summary>
    /// Read QualitySequences from a file
    /// </summary>
     public class QualityFileReader
     {
         StreamReader sr;
         private string name;
         ISequenceValidator validator = new SequenceValidator_AllValid();
         private readonly bool ignoreInvalid = false;
         private readonly bool trimExcess = false;
         ISequenceContainer seqTemplate = SequenceFactory.GetDefaultSequence();


         public QualityFileReader(StreamReader sr)
         {
             this.sr = sr;
         }
         public QualityFileReader(StreamReader sr, ISequenceValidator validator)
         {
             this.sr = sr;
             this.validator = validator;
         }

         public QualityFileReader(string path, ISequenceValidator validator)
             : this(new StreamReader(path, Encoding.ASCII),validator)
         {
         }

         public QualityFileReader(string path)
             : this(new StreamReader(path, Encoding.ASCII))
         {
         }


         public QualityFileReader(StreamReader sr, ISequenceValidator validator, QualityFileReaderOptions options)
         {
             this.sr = sr;
             this.validator = validator;
             if ((options & QualityFileReaderOptions.IgnoreInvalid) == QualityFileReaderOptions.IgnoreInvalid) ignoreInvalid = true;
             if ((options & QualityFileReaderOptions.TrimExcess) == QualityFileReaderOptions.TrimExcess) trimExcess = true;
         }

         public QualityFileReader(string path, ISequenceValidator validator, QualityFileReaderOptions options)
             : this(new StreamReader(path, Encoding.ASCII), validator,options)
         {
         }

         public List<QualitySequence> GetAllSequences()
         {

             List<QualitySequence> sequences = new List<QualitySequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             QualitySequence temp;
             while ((temp = GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} QualityFileReader; All {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, seqTemplate.GetType().ToString()));
#endif
             return sequences;
         }


         public List<QualitySequence> GetSequences(int number)
         {
             List<QualitySequence> sequences = new List<QualitySequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             QualitySequence temp;
             while (sequences.Count < number && (temp = GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} QualityFileReader; A batch of {1} files loaded in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, sw.ElapsedMilliseconds, sequences.Count <= 0 ? 0 : sequences[0].Length, seqTemplate.GetType().ToString()));
#endif
             return sequences;
         }



         public QualitySequence GetNextSequence()
         {
             /*
             >ES6HDNH02EISLO
             28 28 27 27 34 28 28 27 27 28 23 33 25 26 35 28 27 27 28 28 34 27 27 28 26 27 32 25 26 28 27 34 27 27 27 28 28 27 27 27 27 26 28 35 28 34 28 27 37 33 17 27 27 33 25 28 34 27 28 33 25 28 28 34 27 34 27 25 27 27 27 27 27 27 34 27 28 34 27 27 38 34 22 09 28 28 27 27 27 28 34 28 25 27 27 27 28 33 26 28 28 27 28 28 27 28 27 28 33 26 28 28 27 26 28 27 24 28 28 38 34 22 10 27 25 28 28 27 37 33 15 28 27 27 37 33 16 27 27 27 28 27 27 28 27 27 37 33 16 28 27 26 23 26 33 26 25 32 24 27 35 31 12 35 28 27 27 34 28 26 32 24 36 32 13 31 24 21 31 23 27 35 31 11 28 31 26 03 27 28 25 38 34 24 16 09 2 27 27 23 24 24 28 37 33 20 06 21 20 35 30 11 27 28 37 33 16 27 37 33 16 26 33 26 36 32 13 37 33 16 20 27 27 21 37 33 17 20 23 37 32 14 37 33 15 27 27 18 27 27 37 33 16 18 26 34 28 37 33 15 25 20 28 26 18 25 32 24 21
             */
             string line = "";
             QualitySequence ret = null;
             bool invalid = true;
             ISequenceContainer sequence = seqTemplate.NewSequence();

             while (invalid && (line = sr.ReadLine()) != null)
             {
                 line = line.Trim();
                 if (line == "") continue;

                 if (line.StartsWith(">"))
                 {

                     if (name != "" && sequence.Length != 0)
                     {

                         if (validator.IsValid(sequence))
                         {
                             if (trimExcess) sequence.TrimExcess();
                             ret= new QualitySequence(name, sequence);
                             invalid = false;
                         }
                         else if (!ignoreInvalid) throw new InvalidSequenceException(validator,validator.GetInvalid(sequence),name);
                     }

                     //Create new name and sequence
                     name = line.Substring(1, line.Length - 1);
                     sequence = seqTemplate.NewSequence();
                 }
                 else
                 {
                     string[] division = line.Split(new char[] { ' ' });
                     for (int i = 0; i < division.Length; i++)
                     {
                         if (division[i] != "") sequence.Append((byte)Convert.ToInt32(division[i]));
                     }

                 }
             }
             //Also store last file
             if (name != "" && sequence.Length != 0)
             {
                 if (validator.IsValid(sequence))
                 {
                     if (trimExcess) sequence.TrimExcess();
                     ret=new QualitySequence(name, sequence);
                     invalid = false;
                 }
                 else if (!ignoreInvalid) throw new InvalidSequenceException(validator, validator.GetInvalid(sequence), name);
             }

             return ret;
         }


         /// <summary>
         /// Get or set a sequence prototype
         /// </summary>
         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return this.seqTemplate;
             }
             set
             {
                 this.seqTemplate = value;
             }
         }

         public void Close()
         {
             this.sr.Close();
         }

     }

     /// <summary>
     /// Write QualitySequences into a file
     /// </summary>
     public class QualityFileWriter
     {
         private StreamWriter sw;
         int leng = 50;
         public QualityFileWriter(string path) : this(new StreamWriter(path, false, Encoding.ASCII)) { }
         public QualityFileWriter(StreamWriter sw)
         {
             this.sw = sw;
         }

         public void WriteSequences(List<QualitySequence> sequences)
         {
#if DEBUG
             System.Diagnostics.Stopwatch stopw = new System.Diagnostics.Stopwatch();
             stopw.Start();
#endif
             foreach (QualitySequence qs in sequences)
             {
                 WriteSequence(qs);
             }

#if DEBUG
             stopw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} QualityFileWriter; {1} files saved in {2} ms; First file contained {3} quality scores; Sequence used: {4}", DateTime.Now.TimeOfDay, sequences.Count, stopw.ElapsedMilliseconds,
                 sequences.Count <= 0 ? 0 : sequences[0].Length, sequences.Count <= 0 ? "no sequence" : sequences[0].Sequence.GetType().ToString()));


#endif
         }

         public void WriteSequence(QualitySequence sequence)
         {
                 int count = 0;
                 // string line = "";
                 ISequenceContainer seq = sequence.Sequence;
                 sw.WriteLine('>' + sequence.Name);

                 while (seq.Length >= count + leng)
                 {
                     //seq.ToString(count, leng)
                     StringBuilder sb = new StringBuilder("");
                     for (int i = 0; i < leng; i++)
                     {
                         sb.Append((int)seq[count + i]);
                         if (i < leng - 1) sb.Append(' ');
                     }
                     sw.WriteLine(sb.ToString());
                     count += leng;

                 }
                 StringBuilder sb2 = new StringBuilder();
                 for (int i = 0; i < seq.Length - count; i++)
                 {
                     sb2.Append((int)seq[count + i]);
                     if (i < seq.Length - count - 1) sb2.Append(' ');
                 }
                 sw.WriteLine(sb2.ToString());
                 sw.WriteLine();
         }

         public void Close()
         {
             this.sw.Close();
         }


     }


     /// <summary>
     /// Writes a collection of NucSeqInfos into a file
     /// </summary>
     public class NucSeqInfoWriter : INucSeqInfoWriter
     {
         private StreamWriter sw;
         private bool introductionWriten = false;

         public NucSeqInfoWriter(StreamWriter sw)
         {
             this.sw = sw;
         }

         public NucSeqInfoWriter(string path)
             : this(new StreamWriter(path, false, Encoding.ASCII))
         {
         }

         public void WriteSeqInfo(INucSeqInfo sequenceInfo)
         {
             if (!introductionWriten)
             {
                 sw.WriteLine("%>Sequence names");
                 sw.WriteLine("%#Comment");
                 sw.WriteLine("%nfo:\tLength\tIsUpperCase\tParentName\tStart\tEnd\tCount_A\tCount_T\tCount_U\tCount_C\tCount_G\tCount_N\tCount_Special\tCount_InvalidWOIndel\tCount_Indel\tIndividual_sequences");
                 sw.WriteLine();
                 introductionWriten = true;
             }

             for (int i = 0; i < sequenceInfo.Names.Count; i++)
             {
                 sw.WriteLine(String.Format(">{0}", sequenceInfo.Names[i]));
             }

             if (sequenceInfo.Comment != null && sequenceInfo.Comment != "") sw.WriteLine(String.Format("#{0}", sequenceInfo.Comment));

             //                                                                                                                 0        1           2           3           4       5           6           7           8           9           10          11                  12                  13
             sw.WriteLine(String.Format("nfo:\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}", sequenceInfo.Length, sequenceInfo.IsUpper, sequenceInfo.ParentName, sequenceInfo.Start, sequenceInfo.End, sequenceInfo.Count_A, sequenceInfo.Count_T, sequenceInfo.Count_U, sequenceInfo.Count_C, sequenceInfo.Count_G, sequenceInfo.Count_N, sequenceInfo.Count_Special, sequenceInfo.Count_Invalid_woIndel,sequenceInfo.Count_Indel, sequenceInfo.SequenceCount));
             sw.WriteLine();
         }

         public void WriteSeqInfos(List<INucSeqInfo> sequenceInfos)
         {


             foreach (INucSeqInfo ns in sequenceInfos)
             {
                 WriteSeqInfo(ns);
             }
         }

         public void Close()
         {
             this.sw.Close();
         }


     }


    /// <summary>
    /// Write nucleotide sequence infos into a file
    /// </summary>
     public class NucSeqInfoReader : INucSeqInfoReader
     {
         StreamReader sr;
         public NucSeqInfoReader(StreamReader sr)
         {
             this.sr = sr;
         }
         public NucSeqInfoReader(string path)
             : this(new StreamReader(path, Encoding.ASCII))
         {
         }


         public List<INucSeqInfo> GetSequences(int number)
         {
             List<INucSeqInfo> sequences = new List<INucSeqInfo>();

             INucSeqInfo temp;
             while (sequences.Count<number && (temp = GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }
             return sequences;
         }

         public List<INucSeqInfo> GetAllSequences()
         {
             List<INucSeqInfo> sequences = new List<INucSeqInfo>();

             INucSeqInfo temp;
             while ((temp = GetNextSequence()) != null)
             {
                 sequences.Add(temp);
             }

             return sequences;
         }


         public INucSeqInfo GetNextSequence()
         {
             string line = "";
             INucSeqInfo ns = new NucSeqInfo();
             INucSeqInfo ret = null;
             bool invalid = true;

             //bool stored =false;
             while (invalid && (line = sr.ReadLine()) != null)
             {
                 line = line.Trim();
                 if (line == "" || line.StartsWith("%")) continue;

                 if (line.StartsWith(">"))
                 {
                     if (ns.Names == null) ns.Name = line.Substring(1, line.Length - 1);
                     else ns.Names.Add(line.Substring(1, line.Length - 1));
                 }
                 else if (line.StartsWith("#"))
                 {
                     ns.Comment = line.Substring(1, line.Length - 1);
                 }
                 else if (line.StartsWith("nfo:"))
                 {
                     //0     1       2           3           4       5   6       7       8       9       10          11      12              13              14         15                   

                     //nfo:	Length	IsUpperCase	ParentName	Start	End	Count_A	Count_T	Count_U	Count_C	Count_G	Count_N	Count_Special	Count_Invalid	Count_Indel Individual_sequences

                     string[] splits = line.Split(new char[] { '\t' });
                     ns.Length = Convert.ToInt32(splits[1]);
                     ns.IsUpper = splits[2] == "True" ? true : false;
                     ns.ParentName = splits[3];

                     ns.Start = (splits[4] == "" || splits[4] == "null") ? (int?)null : Convert.ToInt32(splits[4]);
                     ns.End = (splits[5] == "" || splits[5] == "null") ? (int?)null : Convert.ToInt32(splits[5]);
                     ns.Count_A = Convert.ToInt32(splits[6]);
                     ns.Count_T = Convert.ToInt32(splits[7]);
                     ns.Count_U = Convert.ToInt32(splits[8]);
                     ns.Count_C = Convert.ToInt32(splits[9]);
                     ns.Count_G = Convert.ToInt32(splits[10]);
                     ns.Count_N = Convert.ToInt32(splits[11]);
                     ns.Count_Special = Convert.ToInt32(splits[12]);
                     ns.Count_Invalid_woIndel = Convert.ToInt32(splits[13]);
                     ns.Count_Indel = Convert.ToInt32(splits[14]);
                     ns.SequenceCount = Convert.ToInt32(splits[15]);

                     ret = ns;
                     invalid = false;
                 }
                 else throw new ArgumentOutOfRangeException("Somethings is wrong with the nucseqinfo");
             }

             //if (!stored)return new List<INucSeqInfo>();
             return ret;
         }

         public void Close()
         {
             this.sr.Close();

         }

     }


     /// <summary>
     /// Decorater for SequenceReader derived classes,
     /// makes sure that each Sequence is uppercase, converts each sequence into uppercase characters
     /// </summary>
     public class NucleotideSequenceReaderDecorator_ObtainQualityScores<TNucleotideSequence> : ISequenceReader<TNucleotideSequence>
         where TNucleotideSequence : INucleotideSequence
     {
         private ISequenceReader<TNucleotideSequence> reader;
         private string[] qualFiles;
         private int bufferQualityFiles;
         private bool ignoreEmptyQuality = false;

         public NucleotideSequenceReaderDecorator_ObtainQualityScores(ISequenceReader<TNucleotideSequence> reader, string[] qualityFiles)
         {

             this.reader = reader;
             this.qualFiles = qualityFiles;
         }

         public NucleotideSequenceReaderDecorator_ObtainQualityScores(ISequenceReader<TNucleotideSequence> reader, string[] qualityFiles, int bufferQualityFiles)
         {
             this.reader = reader;
             this.qualFiles = qualityFiles;
             this.bufferQualityFiles = bufferQualityFiles;
         }

         public NucleotideSequenceReaderDecorator_ObtainQualityScores(ISequenceReader<TNucleotideSequence> reader, string[] qualityFiles, int bufferQualityFiles, bool ignoreEmptyQuality)
         {
             this.reader = reader;
             this.qualFiles = qualityFiles;
             this.bufferQualityFiles = bufferQualityFiles;
             this.ignoreEmptyQuality = ignoreEmptyQuality;
         }


         public List<TNucleotideSequence> GetSequences(int number)
         {
             List<TNucleotideSequence> temp = reader.GetSequences(number);
             return GetQuality(temp);
         }



         public void Close()
         {
             this.reader.Close();
         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return reader.SequencePrototype;
             }
             set
             {
                 reader.SequencePrototype = value;
             }
         }


         public List<TNucleotideSequence> GetAllSequences()
         {
             List<TNucleotideSequence> temp = reader.GetAllSequences();

             return GetQuality(temp);
         }

         private List<TNucleotideSequence> GetQuality(List<TNucleotideSequence> temp)
         {

             Dictionary<string, TNucleotideSequence> nsDict = new Dictionary<string, TNucleotideSequence>(temp.Count);
             foreach (TNucleotideSequence ins in temp)
             {
                 if (!nsDict.ContainsKey(ins.Name))
                 {
                     nsDict.Add(ins.Name, ins);
                 }
                 else throw new InvalidDataException("No two fasta sequences are permited to share the same name (fasta ID)");
             }

             
             int assignedCount = 0;
             for (int i = 0; i < qualFiles.Length; i++)
             {
                 StreamReader sr = new StreamReader(qualFiles[i], Encoding.ASCII);
                 QualityFileReader qr = new QualityFileReader(sr);
                 List<QualitySequence> tempQual;

                 while ((tempQual = qr.GetSequences(bufferQualityFiles)).Count > 0)
                 {
                     foreach (QualitySequence qs in tempQual)
                     {
                         if (nsDict.ContainsKey(qs.Name))
                         {
                             if (nsDict[qs.Name].QualitySequence == null)
                             {
                                 nsDict[qs.Name].QualitySequence = qs;
                                 assignedCount++;
                             }
                             else throw new InvalidOperationException("Nucleotide sequence already contains a quality sequence");

                         }

                         if (assignedCount == nsDict.Count) goto loopEnd;
                     }
                 }
                 tempQual.Clear();
             }
             loopEnd:;

             if (!ignoreEmptyQuality)
             {
                 foreach (TNucleotideSequence seq in temp)
                 {
                     if (seq.QualitySequence == null) throw new ArgumentNullException("Quality sequence of {0} is empty; Make sure that the quality sequence IDs and the fasta IDs are identical", seq.Name);
                 }
             }
             return temp;
         }


         public TNucleotideSequence GetNextSequence()
         {
             throw new Exception("This method is not supported, as it would be too inefficient");
         }

 
     }



     /// <summary>
     /// Creates a number of random sequences, can be plugged into each application instead of an actual fasta reader
     /// </summary>
     public class RandomGenerater_NucleotideSequence : ISequenceReader<NucleotideSequence>
     {
         private int contentGC;
         private int actualCount = 0;
         private int count;
         private int length;
         private ISequenceContainer sequencePrototype = SequenceFactory.GetDefaultSequence();

 


         /// <summary>
         /// Generate a number of sequences with a given CG-content and length, randomly
         /// </summary>
         /// <param name="count">number of sequences to generate</param>
         /// <param name="length">length of the sequences to generate</param>
         /// <param name="contentGC">the GC content of the sequences to generate, number between 1 - 100</param>
         public RandomGenerater_NucleotideSequence(int count, int length, int contentGC)
         {
             this.contentGC = contentGC;
             this.length = length;
             this.count = count;
         }

         /// <summary>
         /// Returns a list of randomly generarted nucleotide sequences
         /// </summary>
         /// <returns></returns>
         public NucleotideSequence GetNextSequence()
         {
             NucleotideSequence ret = null;

                 ISequenceContainer newSeq = SequenceFactory.GetDefaultSequence();

                 Random rand = new Random((int)DateTime.Now.Ticks);
                 //Formula
                 // gcTotal * lengthTotal = gcDone * lengthDone + gcToSet * lengthToSet
                 // => gcToSet = (gcTotal * lengtTotal - gcActual *lengthActual)/lengthToSet
                 double gcActual = 0.0;
                 int gcCount = 0;
                 double gcTotal = (double)contentGC;
                 int gcToSet = (int)gcTotal;


                 for (int i = 0; i < length; i++)
                 {

#if DEBUG
                     int rand1 = rand.Next(1, 101);
                     int rand2 = rand.Next(1, 3);
                     if (rand1 < 1 || rand1 > 100) throw new Exception("fuck random generater");
                     if (rand2 < 1 || rand2 > 2) throw new Exception("they suck");
#endif
                     //lower boundary is inclusive, upper is exclusive
                     if (rand.Next(1, 101) <= gcToSet)
                     {
                         //New G or C
                         if (rand.Next(1, 3) == 1)
                         {
                             newSeq.Append('C');
                         }
                         else
                         {
                             newSeq.Append('G');
                         }

                         //Updated gcActual
                         gcCount++;
                         gcActual = ((double)gcCount) * 100.0 / (double)newSeq.Length;

                     }
                     else
                     {
                         //new A or T
                         if (rand.Next(1, 3) == 1)
                         {
                             newSeq.Append('A');
                         }
                         else
                         {
                             newSeq.Append('T');
                         }

                         //update gcActual
                         gcActual = ((double)gcCount) * 100.0 / (double)newSeq.Length;

                     }

                     gcToSet = (int)(((gcTotal * (double)length) - (gcActual * (double)newSeq.Length)) / ((double)(length - newSeq.Length)));
                 }

                 actualCount++;
                 ret=new NucleotideSequence(String.Format("random_sequence; gc-content_{0}; sequence_number_{1}", contentGC, actualCount), newSeq);


             return ret;

         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return sequencePrototype;
             }
             set
             {
                 this.sequencePrototype = value;
             }
         }

         public void Close()
         {

         }




         public List<NucleotideSequence> GetAllSequences()
         {
             throw new Exception("This method is not supported as the number of sequences to return can be specified in the method GetSequences(number)");
         }



         public List<NucleotideSequence> GetSequences(int number)
         {
             List<NucleotideSequence> ret = new List<NucleotideSequence>();
#if DEBUG
             System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
             sw.Start();
#endif
             while (ret.Count<number)
             {
                 ret.Add(GetNextSequence());
             }

#if DEBUG
             sw.Stop();
             System.Diagnostics.Debug.WriteLine(String.Format("{0} RandomGenerator_NucleotideSequence; {1} sequenceces generated in {2} ms; First file contained {3} nucleotides; Sequence used: {4}", DateTime.Now.TimeOfDay, ret.Count, sw.ElapsedMilliseconds, ret.Count <= 0 ? 0 : ret[0].Length, this.sequencePrototype.GetType().ToString()));

#endif

             return ret;

         }

         
     }


     /// <summary>
     /// Decorater for SequenceReader derived classes,
     /// makes sure that each SimpleSequence has an unique identifier (i.e. Fasta name)
     /// </summary>
     public class SequenceReaderDecorator_UniqueSeqID<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IFastaFileReadable<TSimpleSequence>, new()
     {
         private ISequenceReader<TSimpleSequence> reader;
         private Dictionary<string, int> dict;
         public SequenceReaderDecorator_UniqueSeqID(ISequenceReader<TSimpleSequence> reader)
         {

             dict = new Dictionary<string, int>();
             this.reader = reader;
         }

         public TSimpleSequence GetNextSequence()
         {
             TSimpleSequence seq = reader.GetNextSequence();
             if (seq == null) return default(TSimpleSequence);

                 if (dict.ContainsKey(seq.Name))
                 {
                     string old = seq.Name;
                     int add = dict[seq.Name]++;
                     seq.Name = seq.Name + '_' + add.ToString();

                     //If the new sequence name is per chance already a entry in the dictionary, just add an star until the ID is unique
                     while (dict.ContainsKey(seq.Name))
                     {
                         seq.Name += "*";
                     }
                 }

              
             //The name is unique, add it to the dictionary
             dict.Add(seq.Name, 1);
             return seq;
         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return reader.SequencePrototype;
             }
             set
             {
                 reader.SequencePrototype = value;
             }
         }

         public void Close()
         {
             this.reader.Close();
         }

      

         public List<TSimpleSequence> GetAllSequences()
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;
             while ((temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }
             return toRet;
         }


         public List<TSimpleSequence> GetSequences(int number)
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;

             while (toRet.Count<number && (temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }

             return toRet;
         }
     }


     /// <summary>
     /// Decorator for SequenceReader, extracts for a given collection of IPositionable elements the flanking sequences
     /// May be used in batch mode, reads a sequence and extracts the flanking sequences from the read sequence
     /// Returns an empty collection if no more sequences can be read.
     /// </summary>
     /// <typeparam name="TSimpleSequence"></typeparam>
     public class SequenceReaderDecorator_ExtractFlankingSequences<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IPositionable, IFastaFileReadable<TSimpleSequence>, new()
     {
         private ISequenceReader<TSimpleSequence> reader;
         private bool individualFlankMode;
         private int flankUpstream = 0;
         private int flankDownstream = 0;
         private int totalFlankLength = 0;
         private bool addPositions = true;
         private Queue<TSimpleSequence> toReturn = new Queue<TSimpleSequence>();

         private TSimpleSequence templateSimpleSequence = new TSimpleSequence();

         private ISequenceContainer template;
         private Dictionary<string, List<IPositionable>> toFlank;
         /// <summary>
         /// Extracts flanking sequences for a given collection of elements implementing IPositionable (e.g.: SSRs)
         /// Identification of the sequence in which to search the flanking sequence is done through the fasta ID, so make sure that the .ParentName of the IPositionable fits the Name of the fasta ID. 
         /// Therefore decorators changing the fasta IDs should not be used prior to the extraction process.
         /// The flanking length upstream (5') and downstream (3') of the IPositionable element have to be specified, thus flanking sequences will have variing length. Dependend on the length of the IPositionable element
         /// The total flanking length might be less than the sequence which to flank
         /// </summary>
         public SequenceReaderDecorator_ExtractFlankingSequences(ISequenceReader<TSimpleSequence> reader, List<IPositionable> extractFrom, int flankUpstream, int flankDownstream)
         {
             individualFlankMode = true;
             this.reader = reader;
             toFlank = GroupPositionablesWithEqualParent(extractFrom);
             this.flankDownstream = flankDownstream;
             this.flankUpstream = flankUpstream;
             template = SequenceFactory.GetMediumSequence();
         }

         /// <summary>
         /// Extracts flanking sequences for a given collection of elements implementing IPositionable (e.g.: SSRs)
         /// Identification of the sequence in which to search the flanking sequence is done through the fasta ID, so make sure that the .ParentName of the IPositionable fits the Name of the fasta ID. 
         /// Therefore decorators changing the fasta IDs should not be used prior to the extraction process.
         /// The total length of the flanking sequence has to be specified, thus flanking sequences will always have the same length
         /// The total flanking length might be less than the sequence which to flank
         /// </summary>
         public SequenceReaderDecorator_ExtractFlankingSequences(ISequenceReader<TSimpleSequence> reader, List<IPositionable> extractFrom, int totalFlankLength)
         {
             individualFlankMode = false;
             this.totalFlankLength = totalFlankLength;
             this.reader = reader;
             toFlank = GroupPositionablesWithEqualParent(extractFrom);
             template = SequenceFactory.GetMediumSequence();
         }

         /// <summary>
         /// Extracts flanking sequences for a given collection of elements implementing IPositionable (e.g.: SSRs)
         /// Identification of the sequence in which to search the flanking sequence is done through the fasta ID, so make sure that the .ParentName of the IPositionable fits the Name of the fasta ID. 
         /// Therefore decorators changing the fasta IDs should not be used prior to the extraction process.
         /// The flanking length upstream (5') and downstream (3') of the IPositionable element have to be specified, thus flanking sequences will have variing length. Dependend on the length of the IPositionable element
         /// The total flanking length might be less than the sequence which to flank
         /// </summary>
         public SequenceReaderDecorator_ExtractFlankingSequences(ISequenceReader<TSimpleSequence> reader, List<IPositionable> extractFrom, bool addPositions, int flankUpstream, int flankDownstream)
             : this(reader, extractFrom, flankUpstream, flankDownstream)
         {
             this.addPositions = addPositions;

         }

         /// <summary>
         /// Extracts flanking sequences for a given collection of elements implementing IPositionable (e.g.: SSRs)
         /// Identification of the sequence in which to search the flanking sequence is done through the fasta ID, so make sure that the .ParentName of the IPositionable fits the Name of the fasta ID. 
         /// Therefore decorators changing the fasta IDs should not be used prior to the extraction process.
         /// The total length of the flanking sequence has to be specified, thus flanking sequences will always have the same length
         /// The total flanking length might be less than the sequence which to flank
         /// </summary>
         public SequenceReaderDecorator_ExtractFlankingSequences(ISequenceReader<TSimpleSequence> reader, List<IPositionable> extractFrom, bool addPositions, int totalFlankLength)
             : this(reader, extractFrom, totalFlankLength)
         {
             this.addPositions = addPositions;

         }


         public TSimpleSequence GetNextSequence()
         {

             if (toReturn.Count > 0) return toReturn.Dequeue();

             TSimpleSequence parent;
             do
             {
                 parent = reader.GetNextSequence();
                 if (parent==null) return default(TSimpleSequence);//If no more sequences can be read return the empty element

                     //If the Hash table contains the appropriate key, do:
                     if (toFlank.ContainsKey(parent.Name))
                     {
                         List<IPositionable> childs = toFlank[parent.Name];
                         foreach (IPositionable child in childs)
                         {
                             TSimpleSequence flank;
                             int start;
                             int end;
                             int length;
                             string addToName = "";

                             //Calculate the start position and the length, in dependency of the used flanking mode
                             if (!individualFlankMode)
                             {
                                 length = totalFlankLength;
                                 start = child.Start.Value - ((totalFlankLength - (child.End.Value - child.Start.Value + 1)) / 2) - 1;

                             }
                             else
                             {
                                 length = child.End.Value - child.Start.Value + 1 + flankDownstream + flankUpstream;
                                 start = child.Start.Value - flankUpstream - 1;

                                 //123 456 childstart=6 flank upstream=2
                                 //TTT TTACACACAC
                                 //012 3  start=3;
                             }

                             if (length > parent.Length)
                             {
                                 start = 0;
                                 length = (int)parent.Length;
                             }
                             else
                             {

                                 if (start < 0) start = 0;
                                 if (start + length > parent.Length)
                                 {
                                     //01234567890123 L_14
                                     //TTTTTGTGTGTGTG
                                     start = ((int)parent.Length) - length;
                                 }

                             }

                             end = start + length - 1;


                             int startToFlank_Child = 0;
                             int endToFlank_Child = 0;
                             int featureLength = child.End.Value - child.Start.Value + 1;

                             //123 456 child.start=6 flank upstream=2 start=3
                             //TTT TTACACACAC
                             //--- 0123 startToFlan=2 = child.start-start-1

                             //123456 childstart=1;
                             //ACACACA
                             //0123456 startToflank=0;
                             startToFlank_Child = child.Start.Value - start - 1;

                             //TTTTTAAAAATTTTT
                             //01234567890
                             endToFlank_Child = startToFlank_Child + featureLength - 1;

                             if (addPositions) addToName = String.Format("  Startpos_in_parent={0}  Startpos_here={1}  Length={2}", child.Start.Value, startToFlank_Child + 1, featureLength);


                             flank = templateSimpleSequence.NewInstance(child.ParentName + addToName, child.ParentName, template.NewSequence(parent.Sequence.ToString((int)start, (int)length)), start + 1, end + 1);

                             toReturn.Enqueue(flank); //Push the element into the queue

                         }

                 }
             } while (toReturn.Count < 1); //as long as extracted does not contain sequences continue with the procedure

             return toReturn.Dequeue(); //return one element from the queue
         }

         public void Close()
         {
             this.reader.Close();
         }

         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return template;
             }
             set
             {
                 this.template = value;
             }
         }

         public static Dictionary<string, List<IPositionable>> GroupPositionablesWithEqualParent(List<IPositionable> positionables)
         {
             Dictionary<string, List<IPositionable>> groupedSeq = new Dictionary<string, List<IPositionable>>();
             foreach (IPositionable p in positionables)
             {
                 if (groupedSeq.ContainsKey(p.ParentName))
                 {
                     groupedSeq[p.ParentName].Add(p);
                 }
                 else
                 {
                     List<IPositionable> tempList = new List<IPositionable>();
                     tempList.Add(p);
                     groupedSeq.Add(p.ParentName, tempList);
                 }
             }
             return groupedSeq;
         }



  

         public List<TSimpleSequence> GetAllSequences()
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;
             while ((temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }

             return toRet;
         }


         public List<TSimpleSequence> GetSequences(int number)
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;
             while (toRet.Count<number && (temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }

             return toRet;
         }

      
     }



     /// <summary>
     /// Decorator for SequenceReader, masks the sequence tracts corresponding to a given collection of IPositionable with the specified character
     /// </summary>
     /// <typeparam name="TSimpleSequence"></typeparam>
     public class SequenceReaderDecorator_MaskSequence<TSimpleSequence> : ISequenceReader<TSimpleSequence>
         where TSimpleSequence : ISimpleSequence, IPositionable, IFastaFileReadable<TSimpleSequence>, new()
     {
         private ISequenceReader<TSimpleSequence> reader;
         private TSimpleSequence templateSimpleSequence = new TSimpleSequence();
         private Dictionary<string, List<IPositionable>> toMask;
         private char mask = 'N';

         /// <summary>
         /// the positionables in the collection toMask act as templates for masking the corresponding sequence-tracts in the SimpleSequences passed by the reader
         /// </summary>
         public SequenceReaderDecorator_MaskSequence(ISequenceReader<TSimpleSequence> reader, List<IPositionable> toMask, char mask)
         {

             this.reader = reader;
             this.toMask = SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>.GroupPositionablesWithEqualParent(toMask);
             this.mask = mask;
         }

         public SequenceReaderDecorator_MaskSequence(ISequenceReader<TSimpleSequence> reader, List<IPositionable> toMask)
         {

             this.reader = reader;
             this.toMask = SequenceReaderDecorator_ExtractFlankingSequences<NucleotideSequence>.GroupPositionablesWithEqualParent(toMask);

         }

         public TSimpleSequence GetNextSequence()
         {
             TSimpleSequence sequence = reader.GetNextSequence();

             ISequenceContainer m = sequence.Sequence;
             int countPositionables = 0;
             if (toMask.ContainsKey(sequence.Name))
             {
                 List<IPositionable> positionable = toMask[sequence.Name];
                 countPositionables = positionable.Count;

                 foreach (IPositionable pos in positionable)
                 {
                     for (int i = pos.Start.Value; i <= pos.End.Value; i++)
                     {
                         m[i - 1] = mask;
                     }
                 }
                 string comment = sequence.Comment ?? "";
                 sequence.Comment = comment + String.Format(" {0} Positioables have been masked", countPositionables);
             }


             return sequence;

         }


         public ISequenceContainer SequencePrototype
         {
             get
             {
                 return reader.SequencePrototype;
             }
             set
             {
                 this.reader.SequencePrototype = value;
             }
         }

         public void Close()
         {
             this.reader.Close();
         }



         public List<TSimpleSequence> GetAllSequences()
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;
             while ((temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }
             return toRet;
         }


 
         
         public List<TSimpleSequence> GetSequences(int number)
         {
             List<TSimpleSequence> toRet = new List<TSimpleSequence>();
             TSimpleSequence temp;
             while (toRet.Count<number && (temp = GetNextSequence()) != null)
             {
                 toRet.Add(temp);
             }
             return toRet;
         }

     }




}