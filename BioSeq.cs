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
using System.Text.RegularExpressions;

namespace Bio.Seq
{
    public delegate void BioProgressReporter(object sender, string message);

    #region SequenceValidators
    public interface ISequenceValidator
    {
        /// <summary>
        /// Test whether the sequence is valid
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        bool IsValid(ISequenceContainer sequence);


        /// <summary>
        /// Get a List containing a the positions of the invalid characters
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        List<int> GetInvalid(ISequenceContainer sequence);

        bool IsValid(ISeqInfo sequenceInfo);
        bool IsValid(char c);
    }



    /// <summary>
    /// DNA Sequence Validator: ATCG
    /// Checks if a given Sequence solely consists of AaTtCcGg 
    /// </summary>
    public class DNASequenceValidator_ATCG : ISequenceValidator
    {
        public DNASequenceValidator_ATCG()
        {
        }
        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                        break;
                    default: return false;

                }
            }
            return true;
        }

        /// <summary>
        /// Searches whether the sequence contains invalid characters
        /// Returns the position of the invalid character in a 1 based array;
        /// Returns -1 if no invalid character was found
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();

            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                        break;
                    default: invalid.Add(i + 1);
                        break;

                }
            }
            return invalid;
        }


        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_N > 0 || si.Count_Invalid_woIndel > 0 || si.Count_Special > 0 || si.Count_U > 0 || si.Count_Indel > 0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'C':
                case 'c':
                case 'G':
                case 'g': return true;

                default: return false;
            }
        }
    }




    /// <summary>
    /// DNA Sequence Validator: ATCG-
    /// Checks if a given Sequence solely consists of AaTtCcGg-
    /// </summary>
    public class DNASequenceValidator_ATCGIndel : ISequenceValidator
    {
        public DNASequenceValidator_ATCGIndel()
        {
        }
        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case '-':
                        break;
                    default: return false;

                }
            }
            return true;
        }


        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case '-':
                        break;
                    default: invalid.Add(i + 1);
                        break;
                }
            }
            return invalid;
        }


        public bool IsValid(ISeqInfo sequenceInfo)
        {

            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_N > 0 || si.Count_Invalid_woIndel> 0 || si.Count_Special > 0 || si.Count_U > 0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case '-': return true;

                default: return false;
            }
        }
    }



    /// <summary>
    /// DNA Sequence Validator: ATCGN
    /// Checks if a given Sequence solely consists of AaTtCcGgNn
    /// </summary>
    public class DNASequenceValidator_ATCGN : ISequenceValidator
    {
        public DNASequenceValidator_ATCGN()
        {
        }

        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'n':
                    case 'N':

                        break;
                    default: return false;

                }
            }
            return true;
        }

        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'n':
                    case 'N':

                        break;
                    default: invalid.Add( i+1);
                        break;

                }
            }
            return invalid;
        }

        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_Invalid_woIndel > 0 || si.Count_Special > 0 || si.Count_U > 0 || si.Count_Indel>0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case 'n':
                case 'N': return true;

                default: return false;

            }
        }
    }



    /// <summary>
    /// RNA Sequence Validator: ACGU
    /// Checks if a given Sequence solely consists of AaCcGgUu
    /// </summary>
    public class RNASequenceValidator_ACGU : ISequenceValidator
    {
        public RNASequenceValidator_ACGU()
        {
        }

        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'u':
                    case 'U':

                        break;
                    default: return false;

                }
            }
            return true;
        }

        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'u':
                    case 'U':

                        break;
                    default: invalid.Add(i+1);
                        break;

                }
            }

            return invalid;
        }



        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_N > 0 || si.Count_Invalid_woIndel > 0 || si.Count_Special > 0 || si.Count_T > 0||si.Count_Indel>0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {

                case 'A':
                case 'a':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case 'u':
                case 'U': return true;
                default: return false;


            }
        }
    }



    /// <summary>
    /// RNA Sequence Validator: ACGUN
    /// Checks if a given Sequence solely consists of AaCcGgUuNn
    /// </summary>
    public class RNASequenceValidator_ACGUN : ISequenceValidator
    {
        public RNASequenceValidator_ACGUN()
        {
        }

        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'u':
                    case 'U':
                    case 'n':
                    case 'N':

                        break;
                    default: return false;

                }
            }
            return true;
        }


        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'u':
                    case 'U':
                    case 'n':
                    case 'N':

                        break;
                    default: invalid.Add( i+1);
                        break;

                }
            }
            return invalid;
        }



        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_Invalid_woIndel > 0 || si.Count_Special > 0 || si.Count_T > 0||si.Count_Indel>0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'a':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case 'u':
                case 'U':
                case 'n':
                case 'N': return true;
                default: return false;
            }
        }
    }



    /// <summary>
    /// DNA Sequence Validator: ATCGNYRWSKMBDHV
    /// Checks if a given Sequence solely consists of AaTtCcGgNn
    /// </summary>
    public class DNASequenceValidator_ATCGNYRWSKMBDHV : ISequenceValidator
    {
        public DNASequenceValidator_ATCGNYRWSKMBDHV()
        {
        }

        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'n':
                    case 'N':


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
                    case 'V':
                        break;
                    default: return false;

                }
            }
            return true;
        }

        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'n':
                    case 'N':


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
                    case 'V':
                        break;
                    default: invalid.Add( i+1);
                        break;

                }
            }
            return invalid;
        }



        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_Invalid_woIndel > 0 || si.Count_U > 0||si.Count_Indel>0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'C':
                case 'c':
                case 'G':
                case 'g':
                case 'n':
                case 'N':


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
                case 'V': return true;
                default: return false;
            }
        }
    }



    /// <summary>
    /// DNA Sequence Validator: UppercaseATCGNYRWSKMBDHV

    /// </summary>
    public class DNASequenceValidator_UppercaseATCGNYRWSKMBDHV : ISequenceValidator
    {
        public DNASequenceValidator_UppercaseATCGNYRWSKMBDHV()
        {
        }

        public bool IsValid(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'T':
                    case 'C':
                    case 'G':
                    case 'N':


                    case 'Y':   //Pyrmidin (C or T)  
                    case 'R':   //Purin (G or A)
                    case 'W':   //Weich (A or T)
                    case 'S':   //Stark(C or G)
                    case 'K':   //Ketogruppe (G or T)
                    case 'M':   //Aminogruppe (A or C)
                    case 'B':   //Nicht A (C,G,T)
                    case 'D':   //Nicht C (A,T,G)
                    case 'H':   //Nicht G (A,T,C)
                    case 'V':   //Nicht T (A,C,G)
                        break;
                    default: return false;

                }
            }
            return true;
        }

        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            List<int> invalid = new List<int>();
            for (int i = 0; i < sequence.Length; i++)
            {
                switch (sequence[i])
                {
                    case 'A':
                    case 'T':
                    case 'C':
                    case 'G':
                    case 'N':


                    case 'Y':   //Pyrmidin (C or T)  
                    case 'R':   //Purin (G or A)
                    case 'W':   //Weich (A or T)
                    case 'S':   //Stark(C or G)
                    case 'K':   //Ketogruppe (G or T)
                    case 'M':   //Aminogruppe (A or C)
                    case 'B':   //Nicht A (C,G,T)
                    case 'D':   //Nicht C (A,T,G)
                    case 'H':   //Nicht G (A,T,C)
                    case 'V':   //Nicht T (A,C,G)
                        break;
                    default: invalid.Add( i+1);
                        break;

                }
            }
            return invalid;
        }




        public bool IsValid(ISeqInfo sequenceInfo)
        {
            if (sequenceInfo is INucSeqInfo)
            {
                INucSeqInfo si = sequenceInfo as INucSeqInfo;
                if (si.Count_AllChar < 1) return false;

                if (si.Count_Invalid_woIndel > 0 || si.Count_U > 0||si.Count_Indel>0) return false;
                else return true;

            }
            else return false;
        }
        public bool IsValid(char c)
        {
            switch (c)
            {
                case 'A':
                case 'T':
                case 'C':
                case 'G':
                case 'N':


                case 'Y':   //Pyrmidin (C or T)
                //Purin (G or A)
                case 'R':
                //Weich (A or T)
                case 'W':
                case 'S':   //Stark(C or G)

                //Ketogruppe (G or T)
                case 'K':
                //Aminogruppe (A or C)
                case 'M':
                case 'B':   //Nicht A (C,G,T)

                //Nicht C (A,T,G)
                case 'D':
                //Nicht G (A,T,C)
                case 'H':
                //Nicht T (A,C,G)
                case 'V': return true;
                default: return false;
            }
        }
    }


    /// <summary>
    /// Sequence validator: every character is valid
    /// </summary>
    public class SequenceValidator_AllValid : ISequenceValidator
    {
        public bool IsValid(ISequenceContainer sequence)
        {
            return true;
        }

        public List<int> GetInvalid(ISequenceContainer sequence)
        {
            return new List<int>();
        }


        public bool IsValid(ISeqInfo sequenceInfo)
        {
            return true;
        }

        public bool IsValid(char c)
        {
            return true;
        }
    }


    #endregion


    #region Sequence

    /// <summary>
    /// Factory for creating instances of sequences-containers, either nucleotide or protein sequences
    /// Sequence containers should not be be instantiated directly because the implementation of the containers
    /// might change over time, in which case whole programs must be updated manually. This is not the case when using the Sequence factory
    /// here the default sequences are simply updated and the new sequence containers are automatically incorporated in the classes.
    /// </summary>
    public static class SequenceFactory
    {

        private static ISequenceContainer largeSequencePrototype = new Sequence_ByteArray();
        private static ISequenceContainer defaultSequencePrototype = new Sequence_ByteArray();
        private static ISequenceContainer mediumSizeSequencePrototype = new Sequence_MergedString();
        private static ISequenceContainer smallSizeSequencePrototype = new Sequence_String();


        #region SetPrototypes
        /// <summary>
        /// Set prototype for large sized sequences
        /// prototypes are used as templates for copying, all derived sequences will have the same type
        /// </summary>
        /// <param name="largeSequence">the specified sequence will act as a template for creating large sequences</param>
        public static void SetLargeSequencePrototype(ISequenceContainer largeSequence)
        {
            largeSequencePrototype = largeSequence;
        }

        /// <summary>
        /// Set prototype for default sequences
        /// </summary>
        /// <param name="defaultSequence">the specified sequence will act as a template for creating default sequences</param>
        public static void SetDefaultSequencePrototoype(ISequenceContainer defaultSequence)
        {
            defaultSequencePrototype = defaultSequence;
        }

        /// <summary>
        /// Set prototype for medium sized sequences
        /// </summary>
        /// <param name="mediumSequence">the specified sequence will act as a template for creating default sequences</param>
        public static void SetMediumSequences(ISequenceContainer mediumSequence)
        {
            mediumSizeSequencePrototype = mediumSequence;
        }

        /// <summary>
        /// Set prototype for small sequences
        /// </summary>
        /// <param name="smallSequence">the specified sequence will act as a template for creating small sequences</param>
        public static void SetSmallSequence(ISequenceContainer smallSequence)
        {
            smallSizeSequencePrototype = smallSequence;
        }
        #endregion


        /// <summary>
        /// Should be used for sequences  larger 50.000.000 bp
        /// </summary>
        public static ISequenceContainer GetLargeSequence()
        {
            return largeSequencePrototype.NewSequence();
        }

        /// <summary>
        /// Should be used for sequences  larger 50.000.000 bp
        /// </summary>
        public static ISequenceContainer GetLargeSequence(string initialSequence)
        {
            return largeSequencePrototype.NewSequence(initialSequence);
        }


        /// <summary>
        /// Should be used for sequences  smaller 50.000.000
        /// </summary>
        public static ISequenceContainer GetMediumSequence()
        {
            return mediumSizeSequencePrototype.NewSequence();
        }

        /// <summary>
        /// Should be used for sequences  smaller 50.000.000
        /// </summary>
        public static ISequenceContainer GetMediumSequence(string initialSequence)
        {
            return mediumSizeSequencePrototype.NewSequence(initialSequence);
        }

        /// <summary>
        /// Should be used for obtaining the default sequence container
        /// </summary>
        /// <returns></returns>
        public static ISequenceContainer GetDefaultSequence()
        {
            return defaultSequencePrototype.NewSequence();
        }

        public static ISequenceContainer GetDefaultSmithWatermanSequence()
        {
            return new Sequence_String();
        }

        /// <summary>
        /// Should be used for obtaining the default sequence container
        /// </summary>
        /// <returns></returns>
        public static ISequenceContainer GetDefaultSequence(string initialSequence)
        {
            return defaultSequencePrototype.NewSequence(initialSequence);
        }

        /// <summary>
        /// Should be used for small sequences which are finished built
        /// </summary>
        public static ISequenceContainer GetSmallSequence()
        {
            return smallSizeSequencePrototype.NewSequence();
        }

        /// <summary>
        /// Should be used for small sequences which are finished built
        /// </summary>
        public static ISequenceContainer GetSmallSequence(string initialSequence)
        {
            return smallSizeSequencePrototype.NewSequence(initialSequence);
        }




    }



    /// <summary>
    /// specifies that the element is a feature having a length and a feature name
    /// </summary>
    public interface IFeature
    {
        long Length { get;}
        string FeatureName { get;}
    }


    /// <summary>
    /// Interface which containers for sequences have to implement,
    /// offering basic methods for sequence handling i.e. a fast building process, index access and
    /// converting between the containers.
    /// </summary>
    public interface ISequenceContainer
    {
        /// <summary>
        /// Appending a sequence to the container
        /// </summary>
        bool Append(string toAppend);
        bool Append(char toAppend);
        bool Append(byte toAppend);
        bool Append(byte[] toAppend);
        bool Append(char[] toAppend);
        bool Append(ISequenceContainer toAppend);



        /// <summary>
        /// Compare sequences
        /// </summary>
        bool Equals(ISequenceContainer compare);

        /// <summary>
        /// index access, most importand method; a high access speed is required because in an average bioinformatics
        /// program 1.000.000.000 index accesses are common
        /// </summary>
        char this[int position] { get;set;}

        /// <summary>
        /// Create a new sequence of the same type; 
        /// Can be used for implementing the "Prototype" pattern;
        /// </summary>
        ISequenceContainer NewSequence();

        /// <summary>
        /// Create a new sequence of the same type; 
        /// Can be used for implementing the "Prototype" pattern;
        /// </summary>
        ISequenceContainer NewSequence(string initialSequence);

        /// <summary>
        /// Deletes the content of the container
        /// </summary>
        void Clear();

        /// <summary>
        /// Converts the sequence to another type calling the factory methode of the template sequence
        /// e.g: template.NewSequence(this.ToString());
        /// </summary>
        ISequenceContainer ConvertTo(ISequenceContainer template);

        /// <summary>
        /// Converts the sequence to another type; returns a subsequence of the given sequence
        /// </summary>
        ISequenceContainer ConvertTo(ISequenceContainer template, int startPosition, int length);

        /// <summary>
        /// Returns the length of the container
        /// </summary>
        int Length { get;}

        /// <summary>
        /// Gets or sets the capacity of the container
        /// </summary>
        int Capacity { get; set;}

        /// <summary>
        /// Returns a subsequence of the container, the subsequence has the same type as the container
        /// </summary>
        ISequenceContainer SubSequence(int startPosition, int length);

        /// <summary>
        /// Converts the sequence into a string
        /// </summary>
        string ToString(int startPosition, int length);
        string ToString();

        /// <summary>
        /// Calls the Char.ToUpper() method for each entry in the container;
        /// eg: a->A; t->T
        /// </summary>
        void ToUpper();

        /// <summary>
        /// Resize the sequenze container to the actual number of elements.
        /// If the sequence container is based on a string or stringbuilder, resizing is not necessary.
        /// Furthermore, trimming will only be done if the difference between Capacity and a Occupied space
        /// exceeds 90%
        /// </summary>
        void TrimExcess();

    }






    /// <summary>
    /// Container for sequences; both nucleotide and protein sequences are allowed.
    /// This implementation of a sequence container is based on string and StringBuilder, combining the unique advantages of both; that is the fast speed for building (append) 
    /// the string and the fast index access of the string. Nevertheless indexaccess is only half as fast as of the original string.
    /// One character (nucleotide/AS) requires two bytes.
    /// </summary>
    public class Sequence_MergedString : ISequenceContainer
    {
        //Represents the Sequence
        private StringBuilder seqBuild;
        private bool buildmode = true;
        private string seqStr = "";


        public Sequence_MergedString()
        {
            ToBuildMode();
        }

        public Sequence_MergedString(string initialSequence)
        {
            ToBuildMode();
            seqBuild = new StringBuilder(initialSequence);
        }

        public bool Append(string toAppend)
        {
            if (!buildmode) ToBuildMode();
            seqBuild.Append(toAppend);
            return true;
        }

        public bool Append(ISequenceContainer toAppend)
        {
            if (!buildmode) ToBuildMode();
            seqBuild.Append(toAppend.ToString());
            return true;
        }

        public bool Append(char toAppend)
        {
            if (!buildmode) ToBuildMode();
            seqBuild.Append(toAppend);
            return true;
        }

        public bool Append(char[] toAppend)
        {
            if (!buildmode) ToBuildMode();
            seqBuild.Append(toAppend);
            return true;
        }
        public bool Append(byte toAppend)
        {
            if (!buildmode) ToBuildMode();
            seqBuild.Append((char)toAppend);
            return true;
        }

        public bool Append(byte[] toAppend)
        {
            if (!buildmode) ToBuildMode();
            char[] c = new char[toAppend.Length];
            for (int i = 0; i < toAppend.Length; i++)
            {
                c[i] = (char)toAppend[i];
            }
            seqBuild.Append(c);
            return true;
        }

        /// <summary>
        /// Change the used mode of the container, a string should contain the sequence 
        /// </summary>
        private void ToStringMode()
        {
            this.buildmode = false;
            this.seqStr = seqBuild.ToString();
            this.seqBuild = null;
        }


        /// <summary>
        /// Change the used mode of the container a strinbuilder should contain the sequence
        /// </summary>
        private void ToBuildMode()
        {
            this.buildmode = true;
            this.seqBuild = new StringBuilder(this.seqStr);
            this.seqStr = null;
        }

        /// <summary>
        /// Compares two sequences if they are identical
        /// Fast algorithm, first determines if the length are identical,
        /// only if they are of different length all elements will be compared
        /// </summary>
        public bool Equals(ISequenceContainer compare)
        {
            if (buildmode) ToStringMode();
            if (compare == null) throw new ArgumentNullException("compare");

            if (compare.Length != this.Length)
                return false;
            else
            {

                for (int i = 0; i < compare.Length; i++)
                {
                    if (compare[i] != this[i]) return false;
                }
                return true;
            }


        }

        /// <summary>
        /// Gets or sets the capacity of the sequence container
        /// </summary>
        public int Capacity
        {
            get
            {
                if (buildmode) return seqBuild.Capacity;
                else return seqStr.Length;
            }
            set
            {
                if (buildmode) seqBuild.Capacity = value;
            }
        }

        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation
        /// </summary>
        public ISequenceContainer NewSequence()
        {
            return new Sequence_MergedString();
        }

        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation using a initial sequence
        /// </summary>
        public ISequenceContainer NewSequence(string initialSequence)
        {
            return new Sequence_MergedString(initialSequence);
        }

        /// <summary>
        /// Converts the sequence into another ISequence type, calling the internal factory method of toConvert [.NewSequence()]
        /// </summary>
        public ISequenceContainer ConvertTo(ISequenceContainer toConvert)
        {
            return toConvert.NewSequence(this.ToString());
        }


        public ISequenceContainer ConvertTo(ISequenceContainer template, int startPosition, int length)
        {
            return template.NewSequence(this.ToString(startPosition, length));
        }



        /// <summary>
        /// Gets the actual length of the sequence in the container
        /// </summary>
        public int Length
        {
            get
            {
                if (buildmode) return seqBuild.Length;
                else return seqStr.Length;
            }
        }



        /// <summary>
        /// Returns a subsequence of the sequence, using the same container class [Sequence_MergedString()]
        /// </summary>
        public ISequenceContainer SubSequence(int startPosition, int length)
        {
            if (buildmode) ToStringMode();
            return new Sequence_MergedString(seqStr.Substring(startPosition, length));
        }

        /// <summary>
        /// Returns a string representing a subsequence of the sequence.
        /// </summary>
        public string ToString(int startPosition, int length)
        {
            if (buildmode) ToStringMode();
            return seqStr.Substring(startPosition, length);
        }

        /// <summary>
        /// Return the sequence as a string
        /// </summary>
        public override string ToString()
        {
            if (buildmode) ToStringMode();
            return seqStr;
        }



        public char this[int index]
        {
            get
            {
                if (buildmode) ToStringMode();
                return seqStr[index];
            }
            set
            {
                if (!buildmode) ToBuildMode();
                seqBuild[index] = value;
            }
        }

        /// <summary>
        /// Delets the sequence, and allows creating of a new sequence.
        /// </summary>
        public void Clear()
        {

            seqBuild = new StringBuilder("");
            seqStr = "";
            System.GC.Collect();
            buildmode = true;
        }

        /// <summary>
        /// Converts all characters in the sequence to upper characters.
        /// e.g.: a -> A
        /// </summary>
        public void ToUpper()
        {
            if (!buildmode) this.ToBuildMode();
            for (int i = 0; i < this.Length; i++)
            {
                seqBuild[i] = Char.ToUpperInvariant(seqBuild[i]);
            }
        }




        public void TrimExcess()
        {
            if (buildmode) ToStringMode();
        }

    }

    /// <summary>
    /// Container for sequences; both nucleotide and protein sequences are allowed.
    /// Container is based on a byte array, it allows storing of large sequences (human chromosomes) with a low memory consumption.
    /// Building time is fast, and index access is relatively slow, if a fast index access is required and enough memory is available
    /// the sequence should be converted to the class Sequence_MergedString
    /// One character (nucleotide/AS) requires one bytes.
    /// </summary>
    public class Sequence_ByteArray : ISequenceContainer
    {
        //Represents the Sequence
        private byte[] seq;
        private int position = 0;


        public Sequence_ByteArray()
        {
            seq = new byte[1];
        }

        public Sequence_ByteArray(string initialSequence)
            : this()
        {
            if (initialSequence.Length >= 1) seq = new byte[initialSequence.Length];
            this.Append(initialSequence);
        }

        public bool Append(string toAppend)
        {
            IncreaseCapacity(toAppend.Length);
            for (int i = 0; i < toAppend.Length; i++)
            {
                seq[position] = (byte)toAppend[i];
                position++;
            }
            return true;
        }

        public bool Append(ISequenceContainer toAppend)
        {
            IncreaseCapacity(toAppend.Length);
            for (int i = 0; i < toAppend.Length; i++)
            {
                seq[position] = (byte)toAppend[i];
                position++;
            }
            return true;
        }

        public bool Append(char toAppend)
        {
            IncreaseCapacity(1);
            seq[position] = (byte)toAppend;
            position++;
            return true;
        }
        public bool Append(char[] toAppend)
        {
            IncreaseCapacity(toAppend.Length);
            for (int i = 0; i < toAppend.Length; i++)
            {
                seq[position] = (byte)toAppend[i];
                position++;
            }
            return true;

        }


        public bool Append(byte toAppend)
        {
            IncreaseCapacity(1);
            seq[position] = toAppend;
            position++;
            return true;

        }

        public bool Append(byte[] toAppend)
        {

            IncreaseCapacity(toAppend.Length);
            for (int i = 0; i < toAppend.Length; i++)
            {
                seq[position] = toAppend[i];
                position++;
            }
            return true;
        }

        private void IncreaseCapacity(int toIncrease)
        {
            while (position + 1 + toIncrease > seq.Length)
            {
                byte[] temp = new byte[seq.Length * 2];
                for (int i = 0; i < seq.Length; i++)
                {
                    temp[i] = seq[i];
                }
                seq = temp;
            }
        }
        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation
        /// </summary>
        public ISequenceContainer NewSequence()
        {
            return new Sequence_ByteArray();
        }



        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation using a initial sequence
        /// </summary>
        public ISequenceContainer NewSequence(string initialSequence)
        {
            return new Sequence_ByteArray(initialSequence);
        }




        /// <summary>
        /// Converts the sequence into another ISequence type, calling the internal factory method of toConvert [.NewSequence()]
        /// </summary>
        public ISequenceContainer ConvertTo(ISequenceContainer toConvert)
        {
            return toConvert.NewSequence(this.ToString());
        }

        public ISequenceContainer ConvertTo(ISequenceContainer template, int startPosition, int length)
        {
            //if the template is the same sequence as this container simple call the subsequence method of this container
           // if (template is Sequence_ByteArray) return this.SubSequence(startPosition, length);
            //not necessary as the subsequence method is doing the same thing

            byte[] b = new byte[length];
            ISequenceContainer toRet=template.NewSequence();
            for (int i = 0; i < length; i++)
            {
                b[i] = this.seq[startPosition + i];
            }
            toRet.Append(b);
            return toRet;
        }

        /// <summary>
        /// Compares two sequences if they are identical
        /// Fast algorithm, first determines if the length are identical,
        /// only if they are of different length all elements will be compared
        /// </summary>
        public bool Equals(ISequenceContainer compare)
        {
            if (compare == null) throw new ArgumentNullException("compare");

            if (compare.Length != this.Length)
                return false;
            else
            {

                for (int i = 0; i < compare.Length; i++)
                {
                    if (compare[i] != this[i]) return false;
                }
                return true;
            }


        }

        /// <summary>
        /// Gets or sets the capacity of the sequence container
        /// </summary>
        public int Capacity
        {
            get
            {
                return seq.Length;
            }
            set
            {
                if (value > seq.Length)
                {
                    IncreaseCapacity(value - seq.Length);
                }

            }
        }




        /// <summary>
        /// Gets the actual length of the sequence in the container
        /// </summary>
        public int Length
        {
            get
            {
                return position;
            }
        }


        /// <summary>
        /// Returns a subsequence of the sequence, using the same container class [Sequence_ByteArray()]
        /// </summary>
        public ISequenceContainer SubSequence(int startPosition, int length)
        {
            byte[] b = new byte[length];
            for (int i = 0; i < length; i++)
            {
                b[i] = this.seq[startPosition + i];
            }
            Sequence_ByteArray toRet = new Sequence_ByteArray();
            toRet.Append(b);
            return toRet;
        }

        /// <summary>
        /// Returns a string representing a subsequence of the sequence.
        /// </summary>
        public string ToString(int startPosition, int length)
        {
            StringBuilder sb = new StringBuilder(length);
            for (int i = 0; i < length; i++)
            {
                sb.Append((char)seq[startPosition + i]);
            }

            return sb.ToString();
        }

        /// <summary>
        /// Return the sequence as a string
        /// </summary>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder(position);
            for (int i = 0; i < position; i++)
            {
                sb.Append((char)seq[i]);
            }

            return sb.ToString();
        }



        public char this[int index]
        {
            get
            {

                return (char)seq[index];
            }
            set
            {
                seq[index] = (byte)value;
            }
        }

        /// <summary>
        /// Delets the sequence, and allows creating of a new sequence.
        /// </summary>
        public void Clear()
        {
            seq = new byte[1];
            position = 0;
        }


        /// <summary>
        /// Converts all characters in the sequence to upper characters.
        /// e.g.: a -> A
        /// </summary>
        public void ToUpper()
        {
            for (int i = 0; i < this.Length; i++)
            {
                seq[i] = (byte)Char.ToUpperInvariant((char)seq[i]);
            }

        }




        /// <summary>
        /// Trims the size of the container. Usually the capacity of the container exceeds the length of the sequenze. 
        /// If you are sure that your sequenze is not going to grow anymore and if you want to be memory efficient use this method.
        /// Sequences occupieng more than 90% are of the capacity are per default not trimmed
        /// </summary>
        public void TrimExcess()
        {
            if (this.Capacity * 0.9 > this.Length)
            {
                byte[] temp = new byte[this.Length];
                for (int i = 0; i < this.Length; i++)
                {
                    temp[i] = seq[i];
                }
                seq = temp;
            }
        }


      


  


    }

    /// <summary>
    /// This container should only be used for extremly small sequences which are already finished built. The append methods should not been used.
    /// Container is based solely on the string class. Therefore this class provides a high speed for index based access to its elements,
    /// but has an increadible slow building (append) behaviour. Therefore this class should only contain completed (already built) sequences.
    /// One character (nucleotide/AS) requires two bytes.
    /// </summary>
    public class Sequence_String : ISequenceContainer
    {
        //Represents the Sequence
        private string seq = "";


        public Sequence_String()
        {

        }

        public Sequence_String(string initialSequence)
        {
            seq = initialSequence;
        }

        public bool Append(string toAppend)
        {
            seq += toAppend;
            return true;
        }

        public bool Append(ISequenceContainer toAppend)
        {
            seq += toAppend.ToString();
            return true;
        }

        public bool Append(char toAppend)
        {

            seq += toAppend;
            return true;
        }

        public bool Append(char[] toAppend)
        {
            foreach (char c in toAppend)
            {
                this.Append(c);
            }
            return true;
        }
        public bool Append(byte toAppend)
        {
            seq += (char)toAppend;
            return true;
        }

        public bool Append(byte[] toAppend)
        {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < toAppend.Length; i++)
            {
                sb.Append((char)toAppend[i]);
            }

            seq += sb.ToString();
            return true;
        }

        /// <summary>
        /// Compares two sequences if they are identical
        /// Fast algorithm, first determines if the length are identical,
        /// only if they are of different length all elements will be compared
        /// </summary>
        public bool Equals(ISequenceContainer compare)
        {
            if (compare == null) throw new ArgumentNullException("compare");

            if (compare.Length != this.Length)
                return false;
            else
            {

                for (int i = 0; i < compare.Length; i++)
                {
                    if (compare[i] != this[i]) return false;
                }
                return true;
            }


        }


        /// <summary>
        /// Gets or sets the capacity of the sequence container
        /// </summary>
        public int Capacity
        {
            get
            {
                return seq.Length;
            }
            set
            {

            }
        }



        /// <summary>
        /// Gets the actual length of the sequence in the container
        /// </summary>
        public int Length
        {
            get
            {
                return seq.Length;
            }
        }

        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation
        /// </summary>
        public ISequenceContainer NewSequence()
        {
            return new Sequence_String();
        }

        /// <summary>
        /// Internal factory method for creating an ISequence instance of the same class/implementation with an initial sequence
        /// </summary>
        public ISequenceContainer NewSequence(string initialSequence)
        {
            return new Sequence_String(initialSequence);
        }

        /// <summary>
        /// Returns a subsequence of the sequence, using the same container class [Sequence_ByteArray()]
        /// </summary>
        public ISequenceContainer SubSequence(int startPosition, int length)
        {
            return new Sequence_String(seq.Substring(startPosition, length));
        }

        /// <summary>
        /// Returns a string representing a subsequence of the sequence.
        /// </summary>
        public string ToString(int startPosition, int length)
        {
            return seq.Substring(startPosition, length);
        }

        /// <summary>
        /// Return the sequence as a string
        /// </summary>
        public override string ToString()
        {
            return seq;
        }



        public char this[int index]
        {
            get
            {
                return seq[index];
            }
            set
            {

                //Warning this should not be used
                char[] c = seq.ToCharArray();
                c[index] = value;
                seq = c.ToString();

            }
        }


        /// <summary>
        /// Delets the sequence, and allows creating of a new sequence.
        /// </summary>
        public void Clear()
        {
            seq = "";
        }

        /// <summary>
        /// Converts the sequence into another ISequence type, calling the internal factory method of toConvert [.NewSequence()]
        /// </summary>
        public ISequenceContainer ConvertTo(ISequenceContainer toConvert)
        {
            return toConvert.NewSequence(seq);
        }

        public ISequenceContainer ConvertTo(ISequenceContainer template, int startPosition, int length)
        {
           return template.NewSequence(this.ToString(startPosition, length));
        }


        /// <summary>
        /// Converts all characters in the sequence to upper characters.
        /// e.g.: a -> A
        /// </summary>
        public void ToUpper()
        {
            seq = seq.ToUpperInvariant();
        }



        /// <summary>
        /// Trims excessively occupied space from the sequence container;
        /// This is useless with Sequence_String
        /// </summary>
        public void TrimExcess()
        {
            //should do nothing
        }
    }



    /// <summary>
    /// This container represents a sequence which is mapped to a file in the temporary folder; Therefore the sequence should only occupy a fraction 
    /// of the RAM which is occupied by the common sequences. On the downside the file access is much slower
    /// </summary>
    public class Sequence_File : ISequenceContainer
    {
        //static members
        private static string lastFileName="a";
        private static Regex reg = new Regex(@"(?<first>\w*)(?<second>\w)$", RegexOptions.ExplicitCapture | RegexOptions.Compiled);
        private static Dictionary<string, string> nextExtension;
        private static string tempDirectory;

        //instance members
        private string fullPath;
        private System.IO.FileStream fs;
        private int positionCount = 0;
        private bool positionAtEnd = true;

        /// <summary>
        /// Static constructor will be the first thing called, even before the instance constructor
        /// </summary>
        static Sequence_File()
        {
            tempDirectory=Environment.GetEnvironmentVariable("TEMP");

            nextExtension = new Dictionary<string, string>();
            nextExtension.Add("a", "b");
            nextExtension.Add("b", "c");
            nextExtension.Add("c", "d");
            nextExtension.Add("d", "e");
            nextExtension.Add("e", "f");
            nextExtension.Add("f", "g");
            nextExtension.Add("g", "h");
            nextExtension.Add("h", "i");
            nextExtension.Add("i", "j");
            nextExtension.Add("j", "k");
            nextExtension.Add("k", "l");
            nextExtension.Add("l", "m");
            nextExtension.Add("m", "n");
            nextExtension.Add("n", "o");
            nextExtension.Add("o", "p");
            nextExtension.Add("p", "q");
            nextExtension.Add("q", "r");
            nextExtension.Add("r", "s");
            nextExtension.Add("s", "t");
            nextExtension.Add("t", "u");
            nextExtension.Add("u", "v");
            nextExtension.Add("v", "w");
            nextExtension.Add("w", "x");
            nextExtension.Add("x", "y");
            nextExtension.Add("y", "z");
            nextExtension.Add("z", "aa");

        }


        public Sequence_File()
        {
            Match m=reg.Match(lastFileName);
            //Get new file name
            if (!m.Success) throw new InvalidOperationException("Fatal error in Sequence_File; Regex was not matching");
            string first= m.Groups["first"].Value;
            string second=m.Groups["second"].Value;



            //Create file
            fullPath = System.IO.Path.Combine(tempDirectory,lastFileName );
            fs = new System.IO.FileStream(fullPath, System.IO.FileMode.Create,System.IO.FileAccess.ReadWrite);
            lastFileName = first + nextExtension[second];
        }

        public Sequence_File(string initial_Sequence):this()
        {
            for (int i = 0; i < initial_Sequence.Length; i++)
            {
                fs.WriteByte((byte)initial_Sequence[i]);
            }
        }


        private void SetToEnd()
        {
            fs.Seek(0, System.IO.SeekOrigin.End);
            positionAtEnd = true;
        }
   


        public bool Append(string toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            foreach (char c in toAppend)
            {
                fs.WriteByte((byte)c);
            }
            return true;
        }

        public bool Append(char toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            fs.WriteByte((byte)toAppend);
            return true;
        }



        public bool Append(byte toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            fs.WriteByte(toAppend);
            return true;
        }

        public bool Append(byte[] toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            for (int i = 0; i < toAppend.Length; i++)
            {
                fs.WriteByte(toAppend[i]);
            }
            return true;
        }


        public bool Append(char[] toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            for (int i = 0; i < toAppend.Length; i++)
            {
                fs.WriteByte((byte)toAppend[i]);
            }
            return true;
        }


        public bool Append(ISequenceContainer toAppend)
        {
            if (!positionAtEnd) SetToEnd();
            for (int i = 0; i < toAppend.Length; i++)
            {
                fs.WriteByte((byte)toAppend[i]);
            }
            return true;
        }


        public bool Equals(ISequenceContainer compare)
        {

            if (this.Length != compare.Length) return false;
            for (int i = 0; i < this.Length; i++)
            {
                if (compare[i] != this[i]) return false;
            }
            return true;
        }

        public char this[int position]
        {
            get
            {
                fs.Position=position;
                positionAtEnd = false;
                return (char)fs.ReadByte();

            }
            set
            {
                fs.Position=position;
                fs.WriteByte((byte)value);
                positionAtEnd = false;
            }
        }

        public ISequenceContainer NewSequence()
        {
            return new Sequence_File();
        }

        public ISequenceContainer NewSequence(string initialSequence)
        {
            return new Sequence_File(initialSequence);
        }

        public void Clear()
        {
            fs.Close();
            fs = null;
            fs = new System.IO.FileStream(fullPath, System.IO.FileMode.Create, System.IO.FileAccess.ReadWrite);
            this.positionAtEnd = true;
        }


        public ISequenceContainer ConvertTo(ISequenceContainer template)
        {
            ISequenceContainer cont = template.NewSequence();

            fs.Position = 0;
            byte[] b = new byte[this.Length];
            for(int i=0;i<this.Length;i++)
            {
                b[i] =(byte)fs.ReadByte();
            }
            cont.Append(b);
            return cont;
        }

        public ISequenceContainer ConvertTo(ISequenceContainer template, int startPosition, int length)
        {
            ISequenceContainer cont = template.NewSequence();
            this.positionAtEnd = false;
            fs.Position = startPosition;

            if (template is Sequence_String)
            {
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < length; i++)
                {
                    sb.Append((char)fs.ReadByte());
                }
                cont.Append(sb.ToString());
            }
            else
            {

                byte[] b = new byte[length];
                for (int i = 0; i < length; i++)
                {
                    b[i] = (byte)fs.ReadByte();
                }
                cont.Append(b);
            }
           return cont;
        }


        public int Length
        {
            get { return (int)fs.Length; }
        }

        public int Capacity
        {
            get
            {
                return (int) fs.Length;
            }
            set
            {
            }
        }

        public ISequenceContainer SubSequence(int startPosition, int length)
        {
            byte[] b = new byte[length];
            fs.Position = startPosition;
            for (int i = 0; i < length; i++)
            {
                b[i] = (byte)fs.ReadByte();
            }
            Sequence_File seq = new Sequence_File();
            seq.Append(b);
            return seq;
        }

        public string ToString(int startPosition, int length)
        {
            positionAtEnd = false;

            StringBuilder sb = new StringBuilder();
            fs.Position = startPosition;
            for(int i=0; i<length;i++)
            {
                sb.Append((char)fs.ReadByte());
            }
            return sb.ToString();
        }



        public override string ToString()
        {
            positionAtEnd = false;
            StringBuilder sb = new StringBuilder();
            fs.Position = 0;
            for(int i=0; i<this.Length;i++)
            {
                sb.Append((char)fs.ReadByte());
            }

            return sb.ToString();

        }

        /// <summary>
        /// Converts the sequence into uppercase characters;
        /// Do not call this method of Sequence_File if it is not necessary as this implementation is inefficient
        /// </summary>
        public void ToUpper()
        {
            ISequenceContainer seq = this.ConvertTo(new Sequence_ByteArray());
            seq.ToUpper();
            this.Clear();
            this.Append(seq);
        }

        /// <summary>
        /// Trim excessively used space; this function is useless for Sequence_File implementations
        /// </summary>
        public void TrimExcess()
        {
            
        }

        /// <summary>
        /// Destructor to delete the temporary files
        /// </summary>
        ~Sequence_File()
        {
            fs.Close();
            fs = null;
            System.IO.File.Delete(fullPath);
        }

        /// <summary>
        /// Get the name of the temporary file into which the whole sequence has been stored
        /// </summary>
        public string TemporaryFileName
        {
            get
            {

                return fullPath;
            }
        }
    }
    #endregion


    #region Sequence Interfaces
    public interface ISimpleSequence : IFeature
    {
        string Name { get;set;}
        ISequenceContainer Sequence { get;}
        string Comment { get; set;}
        ISequenceValidator Validator { get; set;}
        bool IsValid();

    }

    public interface IPositionable : IFeature
    {
        int? Start { get;set;}
        int? End { get;set;}
        bool? IsRoot { get;set;}
        string ParentName { get;set;}
    }

    public interface IPositionableSimpleSequence : ISimpleSequence, IPositionable
    {
    }

    public interface INucleotideSequence : IPositionableSimpleSequence
    {
        INucleotideSequence ParentSequence { get;set;}
        INucleotideSequence[] ChildSequence { get;set;}
        QualitySequence QualitySequence { get;set;}
        object Tag { get;set;}

        //Sequence Features
        double GC_Percent { get; }
        int Count_CG { get;}
        int Count_ATU { get;}
        int Count_G { get;}
        int Count_C { get;}
        int Count_A { get;}
        int Count_T { get;}
        //Total number of all N or X
        int Count_N { get;}
        int Count_Special { get;}
        int Count_Indel { get;}
        INucSeqInfo GetNucSeqInfo();

        //etc
    }



    //Provides summary info for a sequence (nucleotide/protein)
    public interface ISeqInfo
    {
        string Name { get;set;}
        List<string> Names { get; set;}
        ISequenceValidator Validator { get; set;}
        bool IsValid();
        string Comment { get; set;}
        /// <summary>
        /// A NucSeqInfo may represent many sequences.
        /// This property indicates the number of sequences upon which this instance is based.
        /// </summary>
        int SequenceCount { get;set;}
    }


    public interface IPosSeqInfo : IPositionable, ISeqInfo
    {
    }

    /// <summary>
    /// Provides a short abbreviation for a nucleotide sequence or may be the summary information for some nucleotide sequences
    /// Does not contain the sequence itself; Contains most of the important information including the summary statistics,
    /// Name, ParentName, End- Start position and so on.
    /// </summary>
    public interface INucSeqInfo : IPosSeqInfo
    {

        bool? IsUpper { get; set;}
        int Count_G { get;set;}
        int Count_C { get;set;}
        int Count_A { get;set;}
        int Count_T { get;set;}
        int Count_U { get; set;}
        int Count_Indel { get;set;}
        //Total number of all N or X
        int Count_N { get;set;}

        /// <summary>
        /// Special characters are: YRWSKMBDHVX
        /// </summary>
        int Count_Special { get;set;}
        int Count_Invalid_woIndel { get; set;}
        long Count_AllChar { get;}

        new long Length { get;set;}             //Hides the inherited Length Property
        int Count_CG { get;}
        int Count_ATU { get;}
        long Count_ATUCG { get;}
        double GC_Percent { get;}

        INucSeqInfo Copy();
        INucSeqInfo Sum(INucSeqInfo toAdd);



    }

    #endregion


    /// <summary>
    /// Base class for each type of nucleotide sequence, independend if it is a local or a global nucleotide sequence
    /// Make sure that the nucleotideSequence variable has a value
    /// </summary>
    public abstract class NucleotideSequence_Base : INucleotideSequence
    {
        private INucleotideSequence parentSequence = null;
        private INucleotideSequence[] childSequence = null;
        private ISequenceContainer nucleotideSequence;
        private int? startPosition = null; //relative to parent
        private int? endPosition = null;   //relative to parent

        private string comment = null;
        private string name = null;
        private string parentName = null;
        private ISequenceValidator validator = null;
        private bool? isValid = null;
        private bool? isRoot = null;
        private object tag = null;

        //Quality sequence
        private QualitySequence qualitySequence;


        public NucleotideSequence_Base()
        {
        }

        public NucleotideSequence_Base(ISequenceContainer nucleotideSequence, string sequenceName)
        {
            this.nucleotideSequence = nucleotideSequence;
            this.name = sequenceName;
        }

        public NucleotideSequence_Base(ISequenceContainer nucleotideSequence, string parentName, int startPosition, int endPosition)
        {
            this.nucleotideSequence = nucleotideSequence;
            this.parentName = parentName;
            this.startPosition = startPosition;
            this.endPosition = endPosition;
        }
        public NucleotideSequence_Base(ISequenceContainer nucleotideSequence, string parentName, string sequenceName, int startPosition, int endPosition)
            : this(nucleotideSequence, sequenceName, startPosition, endPosition)
        {
            this.name = sequenceName;
        }

        public NucleotideSequence_Base(ISequenceContainer nucleotideSequence, INucleotideSequence parentSequence, int startPosition, int endPosition)
            : this(nucleotideSequence, parentSequence.Name, startPosition, endPosition)
        {
            this.parentSequence = parentSequence;
        }
        public NucleotideSequence_Base(ISequenceContainer nucleotideSequence, string sequenceName, INucleotideSequence parentSequence, int startPosition, int endPosition)
            : this(nucleotideSequence, parentSequence, startPosition, endPosition)
        {
            this.name = sequenceName;
        }



        //PROPERTIES
        #region NucleotideSequence_Base Properties

        /// <summary>
        /// The name of the local sequence, many local sequences do not necessarily have a name like microsatellites
        /// </summary>
        public string Name
        {
            get
            {
                return name;
            }
            set
            {
                this.name = value;
            }
        }

        /// <summary>
        /// Comments about the sequence.
        /// Fasta files for instance contain a definition for a comment: the string after the ';'-char
        /// </summary>
        public string Comment
        {
            get
            {
                return this.comment;
            }
            set
            {
                this.comment = value;
            }

        }
        /// <summary>
        /// The name of the global sequence, i.e the sequence of which this sequence is a part of
        /// </summary>
        public string ParentName
        {
            get
            {
                return parentName;
            }
            set
            {
                this.parentName = value;
            }
        }

        /// <summary>
        /// The nucleotide sequence e.g.: AATTCCGGTAGAGACTAGATACCAT
        /// </summary>
        public ISequenceContainer Sequence
        {
            get
            {
                return nucleotideSequence;
            }
            protected set
            {
                this.nucleotideSequence = value;
            }


        }

        /// <summary>
        /// Gets or sets the sequence validator
        /// </summary>
        public ISequenceValidator Validator
        {
            get
            {
                return this.validator;
            }
            set
            {
                this.validator = value;
                this.isValid = null;

            }
        }


        /// <summary>
        /// Determines whether the sequence is valid according the used sequence validator
        /// </summary>
        /// <exception cref="ArgumentNullException">Validator has not been specified</exception>
        public bool IsValid()
        {
            if (validator == null) throw new ArgumentNullException("Could not validate the nucleotide sequence; Validator has not been set");
            if (isValid == null) this.isValid = validator.IsValid(nucleotideSequence);
            return isValid.Value;
        }


        /// <summary>
        /// Is the NucleotideSequence the uppermost parent? The root parent;
        /// The origin of which all other sequences are derived or not?
        /// </summary>
        public bool? IsRoot
        {
            get
            {

                return this.isRoot;
            }
            set
            {
                this.isRoot = value;
            }

        }

        /// <summary>
        /// The global sequence, the sequence of which this sequence is a subsequence
        /// </summary>
        public INucleotideSequence ParentSequence
        {
            get
            {
                return parentSequence;
            }
            set
            {
                this.parentSequence = value;
                if (this.parentSequence != null)
                {
                    this.parentName = parentSequence.Name;
                }
            }
        }

        /// <summary>
        /// The subsequences of this sequence
        /// </summary>
        public INucleotideSequence[] ChildSequence
        {
            get
            {
                return this.childSequence;
            }
            set
            {
                this.childSequence = value;
            }
        }

        /// <summary>
        /// The start position of this sequence compared to the parent sequence, therefore this valua is only valid if a parent sequence exists
        /// The string of the parent sequence is 1-based e.g.:
        /// 12345678 9
        /// AAGGTGAT TTTTTTTTTTT
        /// </summary>
        public virtual int? Start
        {
            get
            {
                return this.startPosition;
            }
            set
            {
                this.startPosition = value;
            }
        }


        /// <summary>
        /// The end position of this sequence compared to the parent sequence, therefore this valua is only valid if a parent sequence exists
        /// The string of the parent sequence is 1-based, and the position of the last nucleotide is denoted e.g.:
        /// 12345678 90123456789->19
        /// AAGGTGAT TTTTTTTTTTT
        /// </summary>
        public virtual int? End
        {
            get
            {
                return this.endPosition;
            }
            set
            {
                this.endPosition = value;
            }

        }

        /// <summary>
        /// The total number of G and C
        /// </summary>
        public int Count_CG
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'G':
                        case 'g':
                        case 'C':
                        case 'c': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of A and T/U
        /// </summary>
        public int Count_ATU
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'A':
                        case 'a':
                        case 'u':
                        case 'U':
                        case 'T':
                        case 't': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of A
        /// </summary>
        public int Count_A
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'A':
                        case 'a': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of T
        /// </summary>
        public int Count_T
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'T':
                        case 't': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of U
        /// </summary>
        public int Count_U
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'U':
                        case 'u': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of C
        /// </summary>
        public int Count_C
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'C':
                        case 'c': count++;
                            break;

                    }

                }
                return count;
            }
        }

        /// <summary>
        /// The total number of G
        /// </summary>
        public int Count_G
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'G':
                        case 'g': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of N
        /// </summary>
        public int Count_N
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'N':
                        case 'n':
                            count++;
                            break;

                    }
                }

                return count;
            }
        }


        /// <summary>
        /// The total number of A
        /// </summary>
        public int Count_Indel
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case '-': count++;
                            break;

                    }
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of Y,R,W,S,K,M,B,D,H or V
        /// </summary>
        public int Count_Special
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'a':
                        case 'A':
                        case 't':
                        case 'T':
                        case 'c':
                        case 'C':
                        case 'g':
                        case 'G':
                        case 'n':
                        case 'N':
                        case 'u':
                        case 'U':
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
                        case 'V': count++;
                            break;


                    }//End Switch
                }
                return count;
            }
        }

        /// <summary>
        /// The total number of invalid characters (characters other than ATCGNU YRWSKMBDHV)
        /// Indels are also considered invalid
        /// </summary>
        public int Count_Invalid
        {
            get
            {
                int count = 0;
                for (int i = 0; i < nucleotideSequence.Length; i++)
                {
                    switch (nucleotideSequence[i])
                    {
                        case 'a':
                        case 'A':
                        case 't':
                        case 'T':
                        case 'c':
                        case 'C':
                        case 'g':
                        case 'G':
                        case 'n':
                        case 'N':
                        case 'u':
                        case 'U':
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
                        case 'V':
                            break;
                        default: count++;
                            break;


                    }//End Switch
                }
                return count;
            }
        }

        /// <summary>
        /// The total length of the nucleotide Sequence
        /// </summary>
        public long Length
        {
            get
            {
                return nucleotideSequence.Length;
            }
        }


        /// <summary>
        /// The GC content of the Sequence in Percent
        /// This value is solely calculated from the number of C or G (not counted are special nucleotides like 'S')
        /// </summary>
        public virtual double GC_Percent
        {
            get
            {
                return 100.0 * ((double)Count_CG / ((double)this.Length-this.Count_Invalid));
            }
        }

        public INucSeqInfo GetNucSeqInfo()
        {
            INucSeqInfo nsi = new NucSeqInfo(this.Name, this.ParentName, this.Comment, this.Start, this.End);

            nsi.IsRoot = this.IsRoot;
            nsi.IsUpper = true;
            nsi.Validator = this.Validator;
            nsi.Length = this.Length;



            for (int i = 0; i < this.nucleotideSequence.Length; i++)
            {
                char c = nucleotideSequence[i];
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
            return nsi;
        }


        /// <summary>
        /// Get or set the quality sequence
        /// </summary>
        public QualitySequence QualitySequence
        {
            get
            {
                return this.qualitySequence;
            }
            set
            {
                if (value.Name != this.Name || value.Length != this.Length) throw new System.IO.InvalidDataException("Not possible to assign quality sequence to the nucleotide sequence; The name and the length of both sequences have to be identical");
                this.qualitySequence = value;
            }
        }

        /// <summary>
        /// Get or set a tag
        /// </summary>
        public object Tag
        {
            get
            {
                return this.tag;
            }
            set
            {
                this.tag = value;
            }
        }

        #endregion




        public virtual string FeatureName
        {
            get { return "Nucleotide Sequence"; }
        }


    }



    /// <summary>
    /// Provides functionality for reading and writing of sequences
    /// for storing sequences and provides information about a sequence
    /// methods for batch renaming and conversion into reverse complement
    /// </summary>
    public class NucleotideSequence : NucleotideSequence_Base, Bio.Seq.IO.IFastaFileReadable<NucleotideSequence>
    {



        private double? gcPercent = null;

        /// <summary>
        /// Empty constructor, nothing will work; 
        /// </summary>
        public NucleotideSequence()
        {

        }

        /// <summary>
        /// Class for a simple nucleotide sequence
        /// per default special nucleotide characters like S,W,M,K are not allowed
        /// </summary>
        public NucleotideSequence(string name, ISequenceContainer nucleotideSequence)
        {
            this.Name = name;
            this.Sequence = nucleotideSequence;

        }

        public NucleotideSequence(string name, string comment, ISequenceContainer nucleotideSequence)
            : this(name, nucleotideSequence)
        {

            this.Comment = comment;

        }

        public NucleotideSequence(string name, string parentName, ISequenceContainer nucleotideSequence, int startPosition, int endPosition)
            : this(name, nucleotideSequence)
        {
            this.ParentName = parentName;
            this.Start = startPosition;
            this.End = endPosition;

        }


        public NucleotideSequence(string name, string parentName, ISequenceContainer nucleotideSequence, int? startPosition, int? endPosition)
            : this(name, nucleotideSequence)
        {
            this.ParentName = parentName;
            this.Start = startPosition;
            this.End = endPosition;

        }

        public NucleotideSequence(string name, string parentName, string comment, ISequenceContainer nucleotideSequence, int startPosition, int endPosition)
            : this(name, parentName, nucleotideSequence, startPosition, endPosition)
        {
            this.Comment = comment;
        }

        public NucleotideSequence(string name, string parentName, string comment, ISequenceContainer nucleotideSequence, int? startPosition, int? endPosition)
            : this(name, parentName, nucleotideSequence, startPosition, endPosition)
        {
            this.Comment = comment;
        }



        /// <summary>
        /// Calculate the gcCount:
        /// +1      for GCS
        /// +0,5    for NXKMYR
        /// +0,66   for VB
        /// +0,33   for HD
        /// </summary>
        private void CalculateGCContent()
        {
            double gcCount = 0.0;

            for (int i = 0; i < this.Sequence.Length; i++)
            {
                switch (Sequence[i])
                {
                    case 'S':   //Stark(C or G)
                    case 's':

                    case 'c':
                    case 'C':
                    case 'g':
                    case 'G': gcCount += 1.0;
                        break;


                    case 'Y':   //Pyrmidin (C or T)
                    case 'y':
                    case 'r':   //Purin (G or A)
                    case 'R':
                    case 'k':   //Ketogruppe (G or T)
                    case 'K':
                    case 'm':   //Aminogruppe (A or C)
                    case 'M':


                    case 'n':
                    case 'N': gcCount += 0.5;
                        break;

                    case 'v':   //Nicht T (A,C,G)
                    case 'V':
                    case 'B':   //Nicht A (C,G,T)
                    case 'b': gcCount += 0.6666666666;
                        break;

                    case 'h':   //Nicht G (A,T,C)
                    case 'H':
                    case 'd':   //Nicht C (A,T,G)
                    case 'D': gcCount += 0.3333333333;
                        break;
                }
            } //end for

            //Calculate GC-content in Percent
            this.gcPercent = ((gcCount * 100.0) / (double)(this.Length-this.Count_Invalid));
        }



        /// <summary>
        /// The GC-content in of the sequence in percent. The value is calculated from the G,C,Y,R,M,X,N,V,B,H,D content
        /// whereas the special characters are weighted according their propability for a G or C
        /// e.g.: 0,3333 for a H (not G)
        /// </summary>
        public override double GC_Percent
        {
            get
            {
                if (this.gcPercent == null)
                {
                    CalculateGCContent();
                }
                return (double)this.gcPercent;
            }
        }

        /// <summary>
        /// Returns a new NucleotideSequence with the given paramters
        /// </summary>
        public NucleotideSequence NewInstance(string name, string comment, ISequenceContainer sequence)
        {

            return new NucleotideSequence(name, comment, sequence);
        }

        //(string name, string parentName, ISequence sequence, int startPosition, int endPosition)
        public NucleotideSequence NewInstance(string name, string parentName, ISequenceContainer sequence, int startPosition, int endPosition)
        {
            return new NucleotideSequence(name, parentName, sequence, startPosition, endPosition);
        }


    }




    #region Exceptions

    /// <summary>
    /// ISequenceValidator: sequence is not valid
    /// ToDo: has still to be tested
    /// </summary>
    public class InvalidSequenceException : ArgumentException
    {
        private ISequenceValidator validator = null;
        private List<int> invalidPositions = new List<int>();
        private string name = "unknown sequence";
        private string message = null;

        public InvalidSequenceException(ISequenceValidator validator, List<int> invalidPositions, string name)
            : base(String.Format("Sequence not valid: {0}", validator.ToString()))
        {
            this.validator = validator;
            this.invalidPositions = invalidPositions;
            this.name = name;
        }
        public InvalidSequenceException(ISequenceValidator validator, List<int> invalidPositions)
            : base(String.Format("Sequence not valid: {0}", validator.ToString()))
        {
            this.validator = validator;
            this.invalidPositions = invalidPositions;
        }
        public InvalidSequenceException(string message)
        {
            this.message = message;
        }

        public ISequenceValidator SequenceValidator
        {
            get
            {
                return this.validator;
            }
        }

        public List<int> InvalidPositions
        {
            get
            {
                return this.invalidPositions;
            }
        }

        public string Name
        {
            get
            {
                return name;
            }
        }
        public string Message
        {
            get
            {
                if (message == null) //if no message is given, create the default message
                {
                    StringBuilder sb = new StringBuilder();
                    sb.AppendLine(String.Format("Invalid sequence as determined by validator: {0}", validator != null ? validator.ToString() : "unknown"));
                    sb.AppendLine(String.Format("Invalid sequence ID: {0}", this.name));
                    sb.AppendLine("Invalid character positions at:");
                    foreach (int i in invalidPositions)
                    {
                        sb.Append(String.Format(" {0}", i));
                    }
                    this.message=sb.ToString();
                }

                return message;
            }
        }

    }

    #endregion


    #region Utility
    /// <summary>
    /// Static methods fore class NucleotideSequence
    /// Provides basic functionality for manipulating NucleotideSequences
    /// </summary>
    public class NucleotideSequenceUtility
    {
        private static Random rand = new Random((int)DateTime.Now.Ticks);


        private static int actualSequencesCount = 0;



        /// <summary>
        /// Method returns reverse complement Fasta-file
        /// if desired with postfix or prefix
        /// </summary>
        /// <param name="fasta">The fasta file to convert</param>
        /// <param name="prePostfix">The Prefix or Postfix</param>
        /// <param name="postfix">Append the postfix variable at the end of the file name?</param>
        /// <returns></returns>
        public static NucleotideSequence GetReverseComplement(NucleotideSequence singleSequence, string prePostfix, bool postfix)
        {
            ISequenceContainer sequence = GetReverseComplement(singleSequence.Sequence);



            //Convert the Name
            string name = singleSequence.Name;
            if (postfix)
            {
                name += prePostfix;
            }
            else
            {
                name = name.Insert(0, prePostfix);
            }

            return new NucleotideSequence(name, sequence);
        }

        /// <summary>
        /// Method returns the reverse complement for a sequence represented through a string
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public static string GetReverseComplement(string sequence)
        {
            StringBuilder trans = new StringBuilder("");
            //Translate Sequence into Reverse Complement
            for (int i = sequence.Length - 1; i >= 0; i--)
            {
                switch (sequence[i])
                {
                    case 'a':
                    case 'A': trans.Append('T');
                        break;
                    case 't':
                    case 'T': trans.Append('A');
                        break;
                    case 'c':
                    case 'C': trans.Append('G');
                        break;
                    case 'g':
                    case 'G': trans.Append('C');
                        break;
                    case 'n':
                    case 'N': trans.Append('N');
                        break;


                    //Special symbols:
                    //the input for this procedure has to be a Sequence object, which proofs if special DNA
                    //symbols are allowed, therefore there is no need to check this

                    case 'Y':   //Pyrmidin (C or T)
                    case 'y': trans.Append('R');
                        break;
                    case 'r':   //Purin (G or A)
                    case 'R': trans.Append('Y');
                        break;
                    case 'w':   //Weich (A or T)
                    case 'W': trans.Append('W');
                        break;
                    case 'S':   //Stark(C or G)
                    case 's': trans.Append('S');
                        break;
                    case 'k':   //Ketogruppe (G or T)
                    case 'K': trans.Append('M');
                        break;
                    case 'm':   //Aminogruppe (A or C)
                    case 'M': trans.Append('K');
                        break;
                    case 'B':   //Nicht A (C,G,T)
                    case 'b': trans.Append('V');
                        break;
                    case 'd':   //Nicht C (A,T,G)
                    case 'D': trans.Append('H');
                        break;
                    case 'h':   //Nicht G (A,T,C)
                    case 'H': trans.Append('D');
                        break;
                    case 'v':   //Nicht T (A,C,G)
                    case 'V': trans.Append('B');
                        break;
                    default: throw new Exception("Sequence contains invalid character, in get reverse complement");

                }
            }
            return trans.ToString();

        }

        /// <summary>
        /// Method returns reverse complement for a sequence represented through an elment implementing the ISequence interface.
        /// </summary>
        /// <param name="sequence">The sequence for which the reverse complement should be computed</param>
        public static ISequenceContainer GetReverseComplement(ISequenceContainer sequence)
        {
            ISequenceContainer trans = sequence.NewSequence();

            //Translate Sequence into Reverse Complement
            for (int i = sequence.Length - 1; i >= 0; i--)
            {
                switch (sequence[i])
                {
                    case 'a':
                    case 'A': trans.Append('T');
                        break;
                    case 't':
                    case 'T': trans.Append('A');
                        break;
                    case 'c':
                    case 'C': trans.Append('G');
                        break;
                    case 'g':
                    case 'G': trans.Append('C');
                        break;
                    case 'x':
                    case 'X':   //Any nucleotide
                    case 'n':
                    case 'N': trans.Append('N');
                        break;


                    //Special symbols:
                    //the input for this procedure has to be a Sequence object, which proofs if special DNA
                    //symbols are allowed, therefore there is no need to check this

                    case 'Y':   //Pyrmidin (C or T)
                    case 'y': trans.Append('R');
                        break;
                    case 'r':   //Purin (G or A)
                    case 'R': trans.Append('Y');
                        break;
                    case 'w':   //Weich (A or T)
                    case 'W': trans.Append('W');
                        break;
                    case 'S':   //Stark(C or G)
                    case 's': trans.Append('S');
                        break;
                    case 'k':   //Ketogruppe (G or T)
                    case 'K': trans.Append('M');
                        break;
                    case 'm':   //Aminogruppe (A or C)
                    case 'M': trans.Append('K');
                        break;
                    case 'B':   //Nicht A (C,G,T)
                    case 'b': trans.Append('V');
                        break;
                    case 'd':   //Nicht C (A,T,G)
                    case 'D': trans.Append('H');
                        break;
                    case 'h':   //Nicht G (A,T,C)
                    case 'H': trans.Append('D');
                        break;
                    case 'v':   //Nicht T (A,C,G)
                    case 'V': trans.Append('B');
                        break;
                    default: throw new Exception("Sequence contains invalid character, in get reverse complement");

                }
            }
            return trans;

        }


        /// <summary>
        /// Method for batch conversion of fasta-files into reverse complement
        /// </summary>
        /// <param name="fastas">The fasta files to transform into the reverse complement</param>
        /// <param name="prePostfix">The Prefix or Postfix</param>
        /// <param name="postfix">Append the postfix variable at the end of the file name?</param>
        /// <returns></returns>
        public static List<NucleotideSequence> GetReverseComplement(List<NucleotideSequence> sequences, string prePostfix, bool postfix)
        {
            List<NucleotideSequence> rcFastas = new List<NucleotideSequence>();
            foreach (NucleotideSequence f in sequences)
            {
                rcFastas.Add(NucleotideSequenceUtility.GetReverseComplement(f, prePostfix, postfix));
            }
            return rcFastas;
        }


        /// <summary>
        /// Renames a collection of fasta-files according to the specified criteria
        /// and returns a new fasta collection
        /// </summary>
        /// <param name="fastas">the fasta collection with the wrong names</param>
        /// <param name="prePostFix">The Post- bzw Prefix for the new names</param>
        /// <param name="postfix">True: append PP-fix at the end of the name</param>
        /// <param name="replaceSpace">True: replace all backspace with underscore</param>
        /// <returns>The renamed fasta collection</returns>
        public static List<NucleotideSequence> RenameSequences(List<NucleotideSequence> sequences, string prePostFix, bool postfix, bool replaceSpace)
        {
            List<NucleotideSequence> renamedFastas = new List<NucleotideSequence>(sequences.Count);
            foreach (NucleotideSequence f in sequences)
            {
                string name = f.Name;
                if (postfix)
                {
                    name += prePostFix;
                }
                else
                {
                    name = name.Insert(0, prePostFix);
                }

                if (replaceSpace)
                {
                    name = name.Replace(' ', '_');
                }
                renamedFastas.Add(new NucleotideSequence(name, f.Sequence));
            }
            return renamedFastas;
        }




        /// <summary>
        /// Returns a randomly created sequence, having the specified properties
        /// </summary>
        /// <param name="lengt"></param>
        /// <param name="gcContent"></param>
        /// <returns></returns>
        public static NucleotideSequence RandomSequence(int length, int gcContent)
        {
            ISequenceContainer newSeq = SequenceFactory.GetDefaultSequence();


            //Formula
            // gcTotal * lengthTotal = gcDone * lengthDone + gcToSet * lengthToSet
            // => gcToSet = (gcTotal * lengtTotal - gcActual *lengthActual)/lengthToSet
            double gcActual = 0.0;
            int gcCount = 0;
            double gcTotal = (double)gcContent;
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
            actualSequencesCount++;
            return new NucleotideSequence(String.Format("Random_sequence_NucleotidSeq;gc-content_{0};sequence_number_{1}", gcContent, actualSequencesCount), newSeq);

        }

        private static char GetMutatedCharacter(char toMutate)
        {
            int character = rand.Next(1, 4);
            switch (toMutate)
            {
                case 'A':
                    switch (character)
                    {
                        case 1: return 'T';
                        case 2: return 'C';
                        case 3: return 'G';
                        default: throw new Exception("this should not happen");

                    }

                case 'T':
                    switch (character)
                    {
                        case 1: return 'A';
                        case 2: return 'C';
                        case 3: return 'G';
                        default: throw new Exception("this should not happen");
                    }

                case 'C':
                    switch (character)
                    {
                        case 1: return 'T';
                        case 2: return 'A';
                        case 3: return 'G';
                        default: throw new Exception("this should not happen");
                    }

                case 'G':
                    switch (character)
                    {
                        case 1: return 'T';
                        case 2: return 'C';
                        case 3: return 'A';
                        default: throw new Exception("this should not happen");
                    }

                default:
                    character = rand.Next(1, 5);
                    switch (character)
                    {
                        case 1: return 'T';
                        case 2: return 'C';
                        case 3: return 'G';
                        case 4: return 'A';
                        default: throw new Exception("this should not happen");
                    }


            }
            // return 'N';
        }

    }


    /// <summary>
    /// Collection of methods for handling ISequence elements
    /// </summary>
    public static class SequenceUtility
    {
        private static Dictionary<string, string> aminoCodons;

        private static Dictionary<string, string> InitializeTranslationHashtable()
        {
            Dictionary<string, string> ac = new Dictionary<string, string>(128);
            ac.Add("UUU", "F");
            ac.Add("UUC", "F");
            ac.Add("UUA", "L");
            ac.Add("UUG", "L");
            ac.Add("UCU", "S");
            ac.Add("UCC", "S");
            ac.Add("UCA", "S");
            ac.Add("UCG", "S");
            ac.Add("UAU", "Y");
            ac.Add("UAC", "Y");
            ac.Add("UAA", "STOP");
            ac.Add("UAG", "STOP");
            ac.Add("UGU", "C");
            ac.Add("UGC", "C");
            ac.Add("UGA", "STOP");
            ac.Add("UGG", "W");

            ac.Add("CUU", "L");
            ac.Add("CUC", "L");
            ac.Add("CUA", "L");
            ac.Add("CUG", "L");
            ac.Add("CCU", "P");
            ac.Add("CCC", "P");
            ac.Add("CCA", "P");
            ac.Add("CCG", "P");
            ac.Add("CAU", "H");
            ac.Add("CAC", "H");
            ac.Add("CAA", "Q");
            ac.Add("CAG", "Q");
            ac.Add("CGU", "R");
            ac.Add("CGC", "R");
            ac.Add("CGA", "R");
            ac.Add("CGG", "R");

            ac.Add("AUU", "I");
            ac.Add("AUC", "I");
            ac.Add("AUA", "I");
            ac.Add("AUG", "M");
            ac.Add("ACU", "T");
            ac.Add("ACC", "T");
            ac.Add("ACA", "T");
            ac.Add("ACG", "T");
            ac.Add("AAU", "N");
            ac.Add("AAC", "N");
            ac.Add("AAA", "K");
            ac.Add("AAG", "K");
            ac.Add("AGU", "S");
            ac.Add("AGC", "S");
            ac.Add("AGA", "R");
            ac.Add("AGG", "R");

            ac.Add("GUU", "V");
            ac.Add("GUC", "V");
            ac.Add("GUA", "V");
            ac.Add("GUG", "V");
            ac.Add("GCU", "A");
            ac.Add("GCC", "A");
            ac.Add("GCA", "A");
            ac.Add("GCG", "A");
            ac.Add("GAU", "D");
            ac.Add("GAC", "D");
            ac.Add("GAA", "E");
            ac.Add("GAG", "E");
            ac.Add("GGU", "G");
            ac.Add("GGC", "G");
            ac.Add("GGA", "G");
            ac.Add("GGG", "G");

            return ac;
        }
        /// <summary>
        /// Indicates whether a specified ISequence consists solely of uppercase letters
        /// </summary>
        public static bool IsUpper(ISequenceContainer sequence)
        {
            for (int i = 0; i < sequence.Length; i++)
            {
                if (Char.IsLower(sequence[i])) return false;
            }
            return true;

        }

        /// <summary>
        /// Reverse the order of the sequence, the function is optimized for Sequence_String and Sequence_ByteArray
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public static ISequenceContainer ReverseSequence(ISequenceContainer sequence)
        {
            ISequenceContainer toRet = sequence.NewSequence();

            //If the sequence is a string use the stringbuilder to construct it
            if (sequence is Sequence_String)
            {
                StringBuilder sb = new StringBuilder();
                sb.Capacity = sequence.Length + 3;
                for (int i = sequence.Length - 1; i >= 0; i--)
                {
                    sb.Append(sequence[i]);
                }
                toRet.Append(sb.ToString());
            }
            else //use a byte array and iterate over the array
            {
                byte[] b = new byte[sequence.Length];
                int k = 0;
                for (int i = sequence.Length - 1; i >= 0; i--)
                {
                    b[k] = (byte)sequence[i];
                    k++;
                }
                toRet.Append(b);
            }

            return toRet;

        }



        /// <summary>
        /// Converts a ISequence element into DNA.
        /// That is U's are replaced with T's
        /// Returns a new instance of ISequence
        /// </summary>
        public static ISequenceContainer ConvertToDNA(ISequenceContainer toConvert)
        {
            ISequenceContainer ret = toConvert.NewSequence();
            ret.Capacity = toConvert.Capacity;

            for (int i = 0; i < toConvert.Length; i++)
            {
                switch (toConvert[i])
                {
                    case 'u': ret.Append('t');
                        break;
                    case 'U': ret.Append('T');
                        break;
                    default: ret.Append(toConvert[i]);
                        break;
                }

            }
            return ret;
        }


        /// <summary>
        /// Converts a ISequence element into RNA
        /// That is T's will be replaced with U's
        /// Returns a new instance of ISequence
        /// </summary>
        public static ISequenceContainer ConvertToRNA(ISequenceContainer toConvert)
        {
            ISequenceContainer ret = toConvert.NewSequence();
            ret.Capacity = toConvert.Capacity;

            for (int i = 0; i < toConvert.Length; i++)
            {
                switch (toConvert[i])
                {
                    case 't': ret.Append('u');
                        break;
                    case 'T': ret.Append('U');
                        break;
                    default: ret.Append(toConvert[i]);
                        break;
                }

            }
            return ret;
        }

        /// <summary>
        /// Converts a ISequence element into a ISequence element of upper characters
        /// Returns a new instance of ISequence.
        /// </summary>
        public static ISequenceContainer ConvertToUpper(ISequenceContainer toConvert)
        {
            ISequenceContainer seq = toConvert.NewSequence();
            for (int i = 0; i < toConvert.Length; i++)
            {
                seq.Append(Char.ToUpperInvariant(toConvert[i]));
            }
            return seq;
        }


        /// <summary>
        /// Converts a RNA/DNA-sequence into a peptide sequence
        /// </summary>
        public static ISequenceContainer ConvertToPeptide(ISequenceContainer toConvert)
        {
            ISequenceContainer ret = SequenceFactory.GetSmallSequence();
            ret.Capacity = toConvert.Capacity / 3;
            toConvert = SequenceUtility.ConvertToRNA(toConvert);
            toConvert.ToUpper();

            if (aminoCodons == null) aminoCodons = InitializeTranslationHashtable();

            bool start = false;
            bool stop = false;
            for (int i = 0; i < toConvert.Length - 2; i += 3)
            {
                string aa = aminoCodons[toConvert.ToString(i, 3)];
                if (aa == "STOP")
                {
                    stop = true;
                    break;
                }
                if (aa == "M") start = true;

                if (start) ret.Append(aa);

            }

            if (!start || !stop) return SequenceFactory.GetSmallSequence("");
            else return ret;

        }


        /// <summary>
        /// All invalid (depends on validator) characters are replaced with a 'N'
        /// The same instance of the sequence is masked
        /// </summary>
        public static void MaskThisSequence(ISequenceContainer toMask, ISequenceValidator validator, char replaceInvalidWith)
        {
            if (!validator.IsValid(replaceInvalidWith)) throw new ArgumentException("Character with which invalid characters should be replaced is itself invalid");
            for (int i = 0; i < toMask.Length; i++)
            {
                if (!validator.IsValid(toMask[i])) toMask[i] = replaceInvalidWith;

            }

        }
        /// <summary>
        /// All invalid (depends on validator) characters are replaced with a 'N'
        /// Returns a new instance of the sequence
        /// </summary>
        /// <exception cref="ArgumentException">Validator has to support 'Nn' as valid characters</exception>
        public static ISequenceContainer MaskSequence(ISequenceContainer toMask, ISequenceValidator validator, char replaceInvalidWith)
        {
            if (!validator.IsValid(replaceInvalidWith)) throw new ArgumentException("Character with which invalid characters should be replaced is itself invalid");

            ISequenceContainer masked = toMask.SubSequence(0, toMask.Length);

            for (int i = 0; i < toMask.Length; i++)
            {
                if (!validator.IsValid(masked[i])) masked[i] = replaceInvalidWith;

            }
            return masked;

        }


        /// <summary>
        /// Creates a hash table for a given ISequence element using the specified step size
        /// </summary>
        public static Dictionary<string, List<int>> GetHashTable(ISequenceContainer toHash, int wordSize)
        {
            //Length:9 stepSize:4 Hashcontent:6
            double maxSize = Math.Pow(4.0, (double)wordSize);
            double steps = ((double)toHash.Length) + 1.0 - ((double)wordSize);
            int size = steps < maxSize ? ((int)(2.0 * steps)) : ((int)(2.0 * maxSize));
            Dictionary<string, List<int>> hash = new Dictionary<string, List<int>>(size);
            int limit = (int)steps;

            for (int i = 0; i < limit; i++)
            {
                string motif = toHash.ToString(i, wordSize);
                if (!hash.ContainsKey(motif))
                {
                    List<int> pos = new List<int>();
                    pos.Add(i + 1);
                    hash.Add(motif, pos);
                }
                else
                {
                    hash[motif].Add(i + 1);
                }
            }
            return hash;
        }


        /// <summary>
        /// Get a DNADictionary for the specified sequence using the given word length
        /// </summary>
        /// <param name="toHash">the sequence for which the DNAdictionary ought to be constructed</param>
        /// <param name="wordSize">the word size</param>
        /// <returns>the DNADictionary containing the positions for each word using a zero based array</returns>
        public static DNADictionary<List<int>> GetDNAHashTable(ISequenceContainer toHash, int wordSize)
        {
            DNADictionary<List<int>> dict = new DNADictionary<List<int>>(wordSize);
            int steps = toHash.Length - wordSize + 1;
            ISequenceContainer working;
            for (int i = 0; i < steps; i++)
            {
                working = toHash.SubSequence(i, wordSize);
                if (dict.ContainsKey(working))
                {
                    dict[working].Add(i);
                }
                else
                {
                    List<int> list = new List<int>();
                    list.Add(i);
                    dict[working] = list;
                }
            }

            return dict;

        }

        /// <summary>
        /// Get a DNADictionary for the specified sequence using the given word length
        /// </summary>
        /// <param name="toHash">the sequence for which the DNAdictionary ought to be constructed</param>
        /// <param name="wordSize">the word size</param>
        /// <returns>the DNADictionary containing the positions for each word using a zero based array</returns>
        public static DNADictionary<List<int>> GetDNAHashTableNonOverlapping(ISequenceContainer toHash, int wordSize)
        {
            DNADictionary<List<int>> dict = new DNADictionary<List<int>>(wordSize);
            int steps = toHash.Length - wordSize + 1;
            ISequenceContainer working;
            for (int i = 0; i < steps; i += wordSize)
            {
                working = toHash.SubSequence(i, wordSize);
                if (dict.ContainsKey(working))
                {
                    dict[working].Add(i);
                }
                else
                {
                    List<int> list = new List<int>();
                    list.Add(i);
                    dict[working] = list;
                }
            }

            return dict;

        }


    }


    


    #endregion


    public class QualitySequence : IPositionableSimpleSequence
    {

        private ISequenceContainer sequence;
        private string name;
        private string parentName;
        private string comment;
        private ISequenceValidator validator;
        private bool? isRoot;
        private int? start=null;
        private int? end=null;

        public QualitySequence(string name, ISequenceContainer sequence)
        {
            this.name = name;
            this.sequence = sequence;
        }

        public QualitySequence(string name, string parentName, ISequenceContainer sequence, int? start, int? end)
        {
            this.name = name;
            this.parentName = parentName;
            this.sequence = sequence;
            this.start = start;
            this.end = end;
        }

        public string Name
        {
            get
            {
                return this.name;
            }
            set
            {
                this.name = value;
            }
        }

        public ISequenceContainer Sequence
        {
            get
            {
                return this.sequence;
            }
        }

        public string Comment
        {
            get
            {
                return this.comment;
            }
            set
            {
                this.comment = value;
            }
        }

        public ISequenceValidator Validator
        {
            get
            {
                return this.validator;
            }
            set
            {
                this.validator = value;
            }
        }

        public bool IsValid()
        {
            return this.validator.IsValid(sequence);
        }


        public long Length
        {
            get
            {
                //  S-----S
                //01234567890 
                return sequence.Length;
            }
        }

        public int? Start
        {
            get
            {
                return this.start;
            }
            set
            {
                this.start = value;
            }
        }

        public int? End
        {
            get
            {
                return this.end;
            }
            set
            {
                this.end = value;
            }
        }

        public bool? IsRoot
        {
            get
            {
                return this.isRoot;
            }
            set
            {
                this.isRoot = value;
            }
        }

        public string ParentName
        {
            get
            {
                return this.parentName;
            }
            set
            {
                this.parentName = value;
            }
        }
        public override string ToString()
        {

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.sequence.Length; i++)
            {
                sb.Append((int)sequence[i]);
                if (i < sequence.Length - 1) sb.Append(' ');
            }

            return sb.ToString();
        }




 


        public string FeatureName
        {
            get { return "Quality Sequence"; }
        }

       

    



    }


    public class NucSeqInfo : INucSeqInfo
    {
        private bool? isUpper = null; //is the sequenc only in upper-case characters.
        private bool? isRoot = null;
        private int? start = null;
        private int? end = null;
        private long length = 0;
        private int a;
        private int t;
        private int u;
        private int g;
        private int c;
        private int n;
        private int indel;
        private int special;
        private int invalid;
        //private int length;
        // private string name="";
        private string parentName = "";
        private string comment = "";
        private ISequenceValidator validator;
        private int fileCount = 1;
        List<string> names = null;

        public NucSeqInfo()
        {
        }
        public NucSeqInfo(string name)
        {
            this.Name = name;
        }
        public NucSeqInfo(string name, string comment)
            : this(name)
        {
            this.comment = comment;
        }
        public NucSeqInfo(string parentName, int startPosition, int endPosition)
        {
            this.parentName = parentName;
            this.start = startPosition;
            this.end = endPosition;

        }
        public NucSeqInfo(string parentName, int? startPosition, int? endPosition)
        {
            this.parentName = parentName;
            this.end = endPosition;
            this.start = startPosition;
        }


        public NucSeqInfo(string name, string parentName, int? startPosition, int? endPosition)
            : this(parentName, startPosition, endPosition)
        {
            this.Name = name;
        }
        public NucSeqInfo(string name, string parentName, int startPosition, int endPosition)
            : this(parentName, startPosition, endPosition)
        {
            this.Name = name;
        }
        public NucSeqInfo(string name, string parentName, string comment, int startPosition, int endPosition)
            : this(name, parentName, startPosition, endPosition)
        {
            this.comment = comment;
        }
        public NucSeqInfo(string name, string parentName, string comment, int? startPosition, int? endPosition)
            : this(name, parentName, startPosition, endPosition)
        {
            this.comment = comment;
        }



        /// <summary>
        /// The name of the sequence
        /// </summary>
        public string Name
        {
            get
            {
                if (this.names != null && names.Count > 0) return names[0];
                return null;
            }
            set
            {
                if (names == null || names.Count == 0)
                {
                    names = new List<string>();
                    names.Add(value);
                }
            }
        }

        /// <summary>
        /// A NucSeqInfo instance may be the sum of many NucSeqInfo, each with its own name.
        /// Therefore if the .FileCount increases >1 the NucSeqInfo will return a list of names.
        /// If the FileCount equals one the property will return the name of the single sequence
        /// </summary>
        public List<string> Names
        {
            get
            {
                if (this.names == null) names = new List<string>();
                return this.names;
            }
            set
            {
                this.names = value;
            }
        }

        /// <summary>
        /// The name of the parent sequence
        /// </summary>
        public string ParentName
        {
            get
            {
                return this.parentName;
            }
            set
            {
                this.parentName = value;
            }
        }


        /// <summary>
        /// Indicates whether the sequence is the root. (Does not have a parent)
        /// </summary>
        public bool? IsRoot
        {
            get
            {
                return this.isRoot;
            }
            set
            {
                this.isRoot = value;
            }
        }

        /// <summary>
        /// The NucSeqInfo may be the summary of many sequences.
        /// This properties indicate the number of sequences on which this
        /// NucSeqInfo instance is based.
        /// If it represents a single sequence this property returns 1.
        /// </summary>
        public int SequenceCount
        {
            get
            {
                return this.fileCount;
            }
            set
            {
                this.fileCount = value;
            }
        }


        /// <summary>
        /// The start position of the sequence with respect to the parent sequence
        /// </summary>
        public int? Start
        {
            get
            {
                return this.start;
            }
            set
            {
                this.start = value;
            }
        }

        /// <summary>
        /// The end position of the sequence with respect to the parent sequence
        /// </summary>
        public int? End
        {
            get
            {
                return this.end;

            }
            set
            {
                this.end = value;
            }
        }


        /// <summary>
        /// The validator used for determining whether the sequence is valid
        /// </summary>
        public ISequenceValidator Validator
        {
            get
            {
                return this.validator;
            }
            set
            {
                this.validator = value;
            }

        }

        /// <summary>
        /// Comments concerning the sequence
        /// </summary>
        public string Comment
        {
            get
            {
                return this.comment;
            }
            set
            {
                this.comment = value;
            }
        }

        /// <summary>
        /// Indicates wether the sequence is valid, using the IValidator specified in the .Validator property
        /// </summary>
        /// <exception cref="ArgumentNullException">NucSeqInfo.IsValid; No validator has been specified</exception>
        public bool IsValid()
        {
            if (validator == null) throw new ArgumentNullException("NucSeqInfo.IsValid; No validator has been specified");
            return validator.IsValid(this);
        }

        /// <summary>
        /// Indicates whether all characters of a sequence are uppercase
        /// </summary>
        public bool? IsUpper
        {
            get
            {
                return this.isUpper;
            }
            set
            {
                this.isUpper = value;
            }
        }
        /// <summary>
        /// The total count of C
        /// </summary>
        public int Count_C
        {
            get
            {
                return this.c;
            }
            set
            {
                this.c = value;
            }
        }

        public int Count_Indel
        {
            get
            {
                return indel;
            }
            set
            {
                this.indel = value; 
            }
        }

        /// <summary>
        /// The total count of G 
        /// </summary>
        public int Count_G
        {
            get
            {
                return this.g;
            }
            set
            {
                this.g = value;
            }
        }

        /// <summary>
        /// The total count of N
        /// </summary>
        public int Count_N
        {
            get
            {
                return this.n;

            }
            set
            {
                this.n = value;
            }
        }

        /// <summary>
        /// The total count of U
        /// </summary>
        public int Count_U
        {
            get
            {
                return this.u;
            }
            set
            {
                this.u = value;
            }
        }

        /// <summary>
        /// The total count of A
        /// </summary>
        public int Count_A
        {
            get
            {
                return this.a;
            }
            set
            {
                this.a = value;
            }
        }

        /// <summary>
        /// The total count of T
        /// </summary>
        public int Count_T
        {
            get
            {
                return this.t;
            }
            set
            {
                this.t = value;
            }
        }

        /// <summary>
        /// The total count of special characters
        /// </summary>
        public int Count_Special
        {
            get
            {
                return this.special;
            }
            set
            {
                this.special = value;
            }
        }

        /// <summary>
        /// The total count of invalid characters; not considering indels
        /// </summary>
        public int Count_Invalid_woIndel
        {
            get
            {
                return this.invalid;
            }
            set
            {
                this.invalid = value;
            }
        }


        /// <summary>
        /// The total count of A,T,U
        /// </summary>
        public int Count_ATU
        {
            get
            {
                return this.Count_A + this.Count_T + this.Count_U;
            }
        }

        /// <summary>
        /// The total count of C,G
        /// </summary>
        public int Count_CG
        {
            get
            {
                return this.Count_C + this.Count_G;
            }
        }

        /// <summary>
        /// The total count of A, T, U, C and G
        /// </summary>
        public long Count_ATUCG
        {
            get
            {
                return ((long)this.Count_ATU) + ((long)this.Count_CG);
            }
        }

        /// <summary>
        /// The total length of the sequence in nt
        /// </summary>
        public long Count_AllChar
        {

            get
            {
                return ((long)this.Count_CG) + ((long)this.Count_Indel)+((long)this.Count_ATU) + ((long)this.Count_Invalid_woIndel) + ((long)this.Count_Special) + ((long)this.Count_N);
            }
        }


        public long Length
        {
            get
            {
                return this.length;
            }
            set
            {
                this.length = value;
            }
        }

        /// <summary>
        /// The GC-content in percent, only G or C are considered (not W or S etc)
        /// </summary>
        public double GC_Percent
        {
            get
            {
                return ((((double)this.Count_CG) * 100.0) / ((double)this.Length));
            }

        }


        /// <summary>
        /// Returns a copy of the NucSeqInfo
        /// </summary>
        public INucSeqInfo Copy()
        {
            INucSeqInfo temp = new NucSeqInfo(this.Name, this.ParentName, this.Comment, this.Start, this.End);
            temp.Count_A = this.Count_A;
            temp.Count_T = this.Count_T;
            temp.Length = this.Length;
            temp.Count_C = this.Count_C;
            temp.Count_G = this.Count_G;
            temp.Count_N = this.Count_N;
            temp.Count_Indel = this.Count_Indel;
            temp.Count_Special = this.Count_Special;
            temp.Count_Invalid_woIndel = this.Count_Invalid_woIndel;
            temp.Count_U = this.Count_U;
            temp.IsRoot = this.IsRoot;
            temp.IsUpper = this.IsUpper;
            temp.Validator = this.Validator;
            temp.SequenceCount = this.SequenceCount;
            temp.Names = this.Names;
            return temp;
        }

        /// <summary>
        /// Create a new NucSeqInfo, which represents the summary of this and the toAdd NucSeqInfo
        /// </summary>
        /// <param name="toAdd">NucSeqInfo to Add to this instance</param>
        /// <returns>the new NucSeqInfo representing the sum of the toAdd and this instance of the NucSeqInfo</returns>
        public INucSeqInfo Sum(INucSeqInfo toAdd)
        {
            INucSeqInfo temp = new NucSeqInfo("");
            temp.Name = null;

            ///Properties which can not be easily summed up
            if (this.ParentName != toAdd.ParentName) temp.ParentName = null;
            else temp.ParentName = this.ParentName;


            List<string> tempNames = new List<string>();
            if (this.Names != null) tempNames.AddRange(this.Names);
            if (toAdd.Names != null) tempNames.AddRange(toAdd.Names);

            temp.Names = tempNames;

            if (this.Comment != toAdd.Comment) temp.Comment = null;
            else temp.Comment = this.Comment;

            if (this.IsRoot != toAdd.IsRoot) temp.IsRoot = null;
            else temp.IsRoot = this.IsRoot;

            if (this.IsUpper != toAdd.IsUpper) temp.IsUpper = null;
            else temp.IsUpper = this.IsUpper;

            if (this.Validator != toAdd.Validator) temp.Validator = null;
            else temp.Validator = this.Validator;


            temp.Start = null;
            temp.End = null;

            ///Properties which be easily summed up
            temp.Length = this.Length + toAdd.Length;
            temp.Count_A = this.Count_A + toAdd.Count_A;
            temp.Count_T = this.Count_T + toAdd.Count_T;
            temp.Count_C = this.Count_C + toAdd.Count_C;
            temp.Count_G = this.Count_G + toAdd.Count_G;
            temp.Count_N = this.Count_N + toAdd.Count_N;
            temp.Count_Indel = this.Count_Indel + toAdd.Count_Indel;
            temp.Count_Special = this.Count_Special + toAdd.Count_Special;
            temp.Count_Invalid_woIndel = this.Count_Invalid_woIndel + toAdd.Count_Invalid_woIndel;
            temp.Count_U = this.Count_U + toAdd.Count_U;
            temp.SequenceCount = this.SequenceCount + toAdd.SequenceCount;

            return temp;

        }


        /// <summary>
        /// Operator for the NucSeqInfo class
        /// </summary>
        /// <param name="first"></param>
        /// <param name="second"></param>
        /// <returns></returns>
        public static NucSeqInfo operator +(NucSeqInfo first, NucSeqInfo second)
        {
            return (NucSeqInfo)first.Sum(second);
        }






        public string FeatureName
        {
            get { return "Nucleotide Sequence"; }
        }



    }

    /// <summary>
    /// Creates a unique value for a given DNA sequence, which has to consist solely of the characters AaTtCcGg
    /// in case the sequence contains another character the hashing algorithm returns null
    /// </summary>
    public class Hash_ATCG
    {
        private int wordSize = 0;
        /// <summary>
        /// Creates a new instance of a Hashing function specific for DNA sequences containing soleley the characters ATCG.
        /// If a sequence contains a character other than these, the algorithm returns null.
        /// A word length has to be specified. (Algorithm returns the same hash value for A, AA, AAA etc)
        /// </summary>
        /// <param name="wordSize">the length of the word for which the hashing class will be used</param>
        public Hash_ATCG(int wordSize)
        {
            if (wordSize > 31) throw new IndexOutOfRangeException("Length of the k-mer is to long, maximum length = 31");
            this.wordSize = wordSize;
        }

        public long? GetHash(string toHash)
        {
            if (toHash.Length != wordSize) throw new ArgumentOutOfRangeException("Length of the word is not equal to the specified length");
            long hash = 0x00;
            for (int i = 0; i < toHash.Length; i++)
            {
                hash = hash << 2;
                switch (toHash[i])
                {
                    case 'A':
                    case 'a': hash = hash | 0x00;
                        break;
                    case 'T':
                    case 't': hash = hash | 0x01;
                        break;
                    case 'c':
                    case 'C': hash = hash | 0x02;
                        break;
                    case 'g':
                    case 'G': hash = hash | 0x03;
                        break;
                    default: return null;
                }

            }
            return hash;

        }

        public long? GetHash(ISequenceContainer toHash)
        {
            if (toHash.Length != wordSize) throw new ArgumentOutOfRangeException("Length of the word is not equal to the specified length");
            long hash = 0x00;
            for (int i = 0; i < toHash.Length; i++)
            {
                hash = hash << 2;
                switch (toHash[i])
                {
                    case 'A':
                    case 'a': hash = hash | 0x00;
                        break;
                    case 'T':
                    case 't': hash = hash | 0x01;
                        break;
                    case 'c':
                    case 'C': hash = hash | 0x02;
                        break;
                    case 'g':
                    case 'G': hash = hash | 0x03;
                        break;
                    default: return null;
                }

            }
            return hash;

        }
    }



    /// <summary>
    /// Represents a lightwight hash table for DNA sequences. 
    /// The hash table will immediatelly allocate an array for all elements which can be stored in the hash table considering four characters and the wordlength (Arraysize= Math.Pow(4,wordlength))
    /// The hash table uses ISequences instances as keys and holds a collection of template T elements.
    /// This DNA hash table uses the Hash_ATCG class for creating the respective hash-values
    /// </summary>
    /// <typeparam name="T"></typeparam>
    public class DNADictionary<T>
    {
        private T[] list;
        private Hash_ATCG hash;
        private int wordsize;
        public DNADictionary(int wordsize)
        {
            if (wordsize > 15) throw new ArgumentOutOfRangeException("the wordsize is restricted to 15, with wordsizes >15 the hash table would allocate too much memory, and the conversion int hash=(long)hash would cause inconsistencies");
            this.wordsize = wordsize;
            list = new T[(int)Math.Pow(4.0, wordsize)];
            hash = new Hash_ATCG(wordsize);
        }


        /// <summary>
        /// Add a value for a given key to the hash table,
        /// if the key is invalid the command will be ignored
        /// </summary>
        /// <param name="key">a DNA sequence which has to consist solely of AaTtCcGg, othewise the sequence is ignored</param>
        /// <param name="toAdd">add this value to the dictionary</param>
        public void Add(ISequenceContainer key, T toAdd)
        {
            int? i = (int?)hash.GetHash(key);
            if (i != null) list[i.Value] = toAdd;
        }

        public void Add(int hash, T toAdd)
        {
            list[hash] = toAdd;
        }

        public T this[ISequenceContainer sequence]
        {
            get
            {
                return list[hash.GetHash(sequence).Value];
            }
            set
            {
                int? i = (int?)hash.GetHash(sequence);
                if (i != null) list[i.Value] = value;
            }
        }

        public T this[int index]
        {
            get
            {
                return list[index];
            }
            set
            {
                list[index] = value;
            }
        }

        public bool ContainsKey(ISequenceContainer key)
        {
            //If the length of the key exceeds the maximal length return false
            if (key.Length != wordsize) return false;
            //Get the hash-value first
            int? h = (int?)hash.GetHash(key);
            //Check if the respective entry in the array is other than null, in that case return true;
            if (h != null && list[h.Value] != null) return true;
            else return false;
        }

        /// <summary>
        /// Returns the total length of the DNADictionary
        /// Attention this value does not represent the actual entries, rather all possible entries, therefore this value is equivalent to 4^wordsize
        /// </summary>
        public int Count
        {
            get
            {
                return list.Length;
            }
        }


        /// <summary>
        /// For a given hashvalue return the corresponding entry in the Dictionary
        /// </summary>
        /// <param name="hash">the hash value, for which the entry ought to be returned</param>
        /// <returns>the entry corresponding to the specified hash value</returns>
        public T GetValueForHash(int hash)
        {
            return list[hash];
        }

        public T[] GetValues()
        {
            return this.list;
        }
    }




}
