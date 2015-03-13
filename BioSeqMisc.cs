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

namespace Bio.Seq.Misc
{
    /// <summary>
    /// Class for assigning unique sequence identifiers for a given collection of sequences (Nucleotide/Protein) which is a prerequisite for
    /// the unambigious identification of sequences and assigning of the appropriate features or traits.
    /// </summary>
    public class CreateUniqueSeqID
    {
        private Dictionary<string, Queue<string>> uniqueIDMatrix;
        public CreateUniqueSeqID(Dictionary<string, Queue<string>> uniqueIDMatrix)
        {
            this.uniqueIDMatrix = uniqueIDMatrix;
        }


        /// <summary>
        /// Constructs a matrix (hash table) for a given list of all sequence identifies. Sequences with identical IDs will be 
        /// grouped together and subsequently unique IDs will be assigned.
        /// A unique ID for a given not unique sequence Identifier migth be retrieved from the hash table where the not unique ID acts as key
        /// and the Queue holds the unique identifiers. Access through the .Dequeue() method
        /// </summary>
        public static Dictionary<string,System.Collections.Generic.Queue<string>> GetUniqueIDMatrix(List<string> sequenceIDs)
        {
            Dictionary<string, List<string>> tempmatrix = new Dictionary<string, List<string>>();

            //First group all identical ids together in a Queue
            for (int i = 0; i < sequenceIDs.Count; i++)
            {
                List<string> activeList;
                if (tempmatrix.ContainsKey(sequenceIDs[i]))
                {
                    activeList = tempmatrix[sequenceIDs[i]];
                    activeList.Add(sequenceIDs[i]);
                }
                else
                {

                    activeList = new List<string>();
                    activeList.Add(sequenceIDs[i]);
                    tempmatrix.Add(sequenceIDs[i], activeList);

                }



            }

            Dictionary<string, Queue<string>> matrix = new Dictionary<string, Queue<string>>(tempmatrix.Count);



            foreach (string key in tempmatrix.Keys)
            {
                List<string> notUnique = tempmatrix[key];
                if (notUnique.Count > 1)
                {
                    for (int i = 0; i < notUnique.Count; i++)
                    {

                        notUnique[i] = notUnique[i] + '_' + (i + 1).ToString();

                        //If per chance the new name is already an entry in the matrix add a '*' for each found identity
                        while (tempmatrix.ContainsKey(notUnique[i]))
                        {
                            notUnique[i] += '*';
                        }
                    }
                }

                Queue<string> uniqueID = new Queue<string>(notUnique.Count);

                foreach (string name in notUnique)
                {
                    uniqueID.Enqueue(name);
                }
                matrix.Add(key, uniqueID);

            }

            return matrix;

        }

        public static Dictionary<string,System.Collections.Generic.Queue<string>> GetUniqueIDMatrix(List<INucleotideSequence> nucleotideSequences)
        {
            List<string> names = new List<string>(nucleotideSequences.Count);
            foreach (NucleotideSequence ns in nucleotideSequences)
            {
                names.Add(ns.Name);
            }
            return GetUniqueIDMatrix(names);

        }

        public static Dictionary<string,System.Collections.Generic.Queue<string>> GetUniqueIDMatrix(List<INucSeqInfo> nucSeqInfos)
        {
            List<string> names = new List<string>(nucSeqInfos.Count);

            foreach (NucSeqInfo nsi in nucSeqInfos)
            {
                names.Add(nsi.Name);

            }

            return GetUniqueIDMatrix(names);
        }

        /// <summary>
        /// Checks if a given collection of sequence identifiers contains only unique identifiers, which is a crucial prerequisite for many applications.
        /// </summary>
        public static bool HasUniqueIDs(List<string> sequenceIDs)
        {
            Dictionary<string, bool> matrix = new Dictionary<string, bool>();
            for (int i = 0; i < sequenceIDs.Count; i++)
            {
                if (matrix.ContainsKey(sequenceIDs[i])) return false;
                else matrix.Add(sequenceIDs[i], true);

            }
            return true;
        }


        /// <summary>
        /// Get an unique identifier for a given sequence name. If ID of a sequence is not unique the program will append _1, _2, _3 etc to the respective sequence names.
        /// The seqID will be used as key in a hash-table which contains the unique IDs (uniqueIDmatrix)
        /// </summary>
        public string GetUniqueID(string seqID)
        {
            if (uniqueIDMatrix.ContainsKey(seqID) && uniqueIDMatrix[seqID].Count > 0)
            {
                return uniqueIDMatrix[seqID].Dequeue();

            }
            else throw new ArgumentException("The sequence matrix does not contain this sequence, make sure to use the appropriate sequence matrix");
        }
    }





    public class NucleotideSequenceMutator
    {
        private struct PositionNucleotidePair
        {
            public PositionNucleotidePair(int position, char nucleotide)
            {
                this.position = position;
                this.nucleotide = nucleotide;
            }
            public int position;
            public char nucleotide;
        }

        private Random rand = new Random((int)DateTime.Now.Ticks);
        private char[] allowed = new char[] { 'A', 'T', 'C', 'G' };

        public NucleotideSequenceMutator()
        {
        }




        /// <summary>
        /// Returns a new instance of a NucleotideSequence, which has been mutated with base substitution.
        /// The given NucleotideSequence acts as a template for the mutation process
        /// </summary>
        /// <param name="mutationTemplate">the template NucleotideSequence, the new instance will resemble this template to the specified degree</param>
        /// <param name="similarity">specifies the degree of similarity of the new NucleotideSequence to the template sequence. value between 1-100 (%)</param>
        /// <returns></returns>
        public INucleotideSequence Mutate_BaseSubstitutions_Similarity(INucleotideSequence mutationTemplate, int similarity)
        {
            int mutations = (int)Math.Round((((100.0 - (double)similarity) / 100.0) * (double)mutationTemplate.Sequence.Length), 0);
            INucleotideSequence ns = this.MutateCount(mutationTemplate, mutations, new char[] { 'A', 'T', 'C', 'G' });
            ns.Name += String.Format("_Similarity {0}%", similarity);
            return ns;

        }



        public INucleotideSequence Mutate_BaseSubstitutions(INucleotideSequence mutationTemplate, int mutationsCount)
        {

            INucleotideSequence ns = this.MutateCount(mutationTemplate, mutationsCount, new char[] { 'A', 'T', 'C', 'G' });
            ns.Name += String.Format("_Mutations {0}", mutationsCount);
            return ns;

        }

        public INucleotideSequence Mutate_Transitions(INucleotideSequence mutationTemplate, int mutationsCount)
        {

            INucleotideSequence ns = this.TransitionCount(mutationTemplate, mutationsCount);
            ns.Name += String.Format("_Mutations {0}", mutationsCount);
            return ns;

        }

        public INucleotideSequence Mutate_Transversions(INucleotideSequence mutationTemplate, int mutationsCount)
        {

            INucleotideSequence ns = this.TransversionCount(mutationTemplate, mutationsCount);
            ns.Name += String.Format("_Mutations {0}", mutationsCount);
            return ns;

        }


        /// <summary>
        /// Returns a new instance of a NucleotideSequence, which has been mutated with base substitution.
        /// The given NucleotideSequence acts as a template for the mutation process
        /// </summary>
        /// <param name="mutationTemplate">the template NucleotideSequence, the new instance will resemble this template to the specified degree</param>
        /// <param name="similarity">specifies the degree of similarity of the new NucleotideSequence to the template sequence. value between 1-100 (%)</param>
        /// <returns></returns>
        public INucleotideSequence Mutate_Indels(INucleotideSequence mutationTemplate, int numberOfIndelEvents, int averageIndelSize)
        {
            //First represent the NucleotidSequence in a new way, a linked list of keyValuePairs
            List<int> positions = new List<int>(mutationTemplate.Sequence.Length);
            Dictionary<int, bool> indelPositions = new Dictionary<int, bool>(numberOfIndelEvents);

            for (int i = 0; i < mutationTemplate.Sequence.Length; i++)
            {
                positions.Add(i);
            }


            //Create the mutations; use a list to store the newly mutated nucleotides

            Bio.Stat.RandomGaussian randG = new Bio.Stat.RandomGaussian();
            for (int i = 0; i < numberOfIndelEvents; i++)
            {
                int posIn = rand.Next(0, positions.Count);
                indelPositions.Add(positions[posIn], true);
                positions.RemoveAt(posIn);
            }


            //Create a new ISequence
            ISequenceContainer retSeq = SequenceFactory.GetDefaultSequence();
            for (int i = 0; i < mutationTemplate.Sequence.Length; i++)
            {
                if (indelPositions.ContainsKey(i))
                {
                    int indelSize = 0;
                    while (((indelSize = (int)randG.NextGaussian((double)averageIndelSize, Math.Sqrt((double)averageIndelSize))) < 1))
                    {
                    }
                    if (rand.Next(1, 3) == 1)
                    {
                        //Deletion
                        i += indelSize;
                    }
                    else
                    {
                        for (int k = 0; k < indelSize; k++)
                        {
                            switch (rand.Next(1, 5))
                            {
                                case 1: retSeq.Append('A');
                                    break;
                                case 2: retSeq.Append('T');
                                    break;
                                case 3: retSeq.Append('C');
                                    break;
                                case 4: retSeq.Append('G');
                                    break;
                                default: throw new Exception("This should not happen");
                            }
                        }
                    }

                }
                else
                {
                    retSeq.Append(mutationTemplate.Sequence[i]);
                }

            }

            return new NucleotideSequence(mutationTemplate.Name + String.Format("_Mutated;_Indels:{0}_Average Size:{1}", numberOfIndelEvents, averageIndelSize), retSeq);

        }

        #region Private Methods
        private char GetTransition(char toMutate)
        {

            switch (toMutate)
            {
                case 'A': return 'G';
                case 'G': return 'A';
                case 'C': return 'T';
                case 'T': return 'C';
                default: throw new Exception();
            }
            // return '\0';
        }

        private char GetTransversion(char toMutate)
        {
            char[] allowed = new char[1];
            switch (toMutate)
            {
                case 'A': allowed = new char[] { 'T', 'C' };
                    break;
                case 'T': allowed = new char[] { 'A', 'G' };
                    break;
                case 'C': allowed = new char[] { 'A', 'G' };
                    break;
                case 'G': allowed = new char[] { 'T', 'G' };
                    break;
                    throw new Exception();
            }
            return allowed[rand.Next(0, 2)];
        }

        private char GetMutatedChar(char toMutate, char[] allowedChars)
        {
            List<char> choice = new List<char>(3);
            foreach (char c in allowedChars)
            {
                if (c != toMutate)
                {
                    choice.Add(c);
                }
            }
            return choice[rand.Next(0, choice.Count)];
        }

        private INucleotideSequence MutateCount(INucleotideSequence mutationTemplate, int mutationCount, char[] allowedChars)
        {
            //First represent the NucleotidSequence in a new way, a linked list of keyValuePairs
            List<PositionNucleotidePair> characters = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);

            for (int i = 0; i < mutationTemplate.Sequence.Length; i++)
            {
                characters.Add(new PositionNucleotidePair(i, mutationTemplate.Sequence[i]));
            }
            //Create the mutations; use a list to store the newly mutated nucleotides

            List<PositionNucleotidePair> newSequence = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);
            //int mutations = (int)Math.Round((((100.0 - (double)similarity) / 100.0) * (double)characters.Count), 0);
            for (int i = 0; i < mutationCount; i++)
            {
                int positionOfMutation = rand.Next(0, characters.Count);
                PositionNucleotidePair toMutate = characters[positionOfMutation];
                characters.RemoveAt(positionOfMutation);

                toMutate.nucleotide = GetMutatedChar(toMutate.nucleotide, allowedChars);
                newSequence.Add(toMutate);
            }


            //Write all other characters into the list containing the new sequence
            for (int i = 0; i < characters.Count; i++)
            {
                newSequence.Add(characters[i]);
            }



            //Sorter using anonymous delegate to sort 
            newSequence.Sort(delegate(PositionNucleotidePair p1, PositionNucleotidePair p2)
        {
            if (p1.position > p2.position) return 1;
            else if (p2.position > p1.position) return -1;
            else return 0;
        }
                );

            //Create a new ISequence
            ISequenceContainer retSeq = SequenceFactory.GetDefaultSequence();
            for (int i = 0; i < newSequence.Count; i++)
            {
                retSeq.Append(newSequence[i].nucleotide);
            }

            return new NucleotideSequence(mutationTemplate.Name + "_Mutated;", retSeq);
        }

        private INucleotideSequence TransitionCount(INucleotideSequence mutationTemplate, int mutationCount)
        {
            //First represent the NucleotidSequence in a new way, a linked list of keyValuePairs
            List<PositionNucleotidePair> characters = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);

            for (int i = 0; i < mutationTemplate.Sequence.Length; i++)
            {
                characters.Add(new PositionNucleotidePair(i, mutationTemplate.Sequence[i]));
            }
            //Create the mutations; use a list to store the newly mutated nucleotides

            List<PositionNucleotidePair> newSequence = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);
            //int mutations = (int)Math.Round((((100.0 - (double)similarity) / 100.0) * (double)characters.Count), 0);
            for (int i = 0; i < mutationCount; i++)
            {
                int positionOfMutation = rand.Next(0, characters.Count);
                PositionNucleotidePair toMutate = characters[positionOfMutation];
                characters.RemoveAt(positionOfMutation);

                toMutate.nucleotide = GetTransition(toMutate.nucleotide);
                newSequence.Add(toMutate);
            }


            //Write all other characters into the list containing the new sequence
            for (int i = 0; i < characters.Count; i++)
            {
                newSequence.Add(characters[i]);
            }



            //Sorter using anonymous delegate to sort 
            newSequence.Sort(delegate(PositionNucleotidePair p1, PositionNucleotidePair p2)
        {
            if (p1.position > p2.position) return 1;
            else if (p2.position > p1.position) return -1;
            else return 0;
        }
                );

            //Create a new ISequence
            ISequenceContainer retSeq = SequenceFactory.GetDefaultSequence();
            for (int i = 0; i < newSequence.Count; i++)
            {
                retSeq.Append(newSequence[i].nucleotide);
            }

            return new NucleotideSequence(mutationTemplate.Name + "_Mutated_Transitions;", retSeq);
        }



        private INucleotideSequence TransversionCount(INucleotideSequence mutationTemplate, int mutationCount)
        {
            //First represent the NucleotidSequence in a new way, a linked list of keyValuePairs
            List<PositionNucleotidePair> characters = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);

            for (int i = 0; i < mutationTemplate.Sequence.Length; i++)
            {
                characters.Add(new PositionNucleotidePair(i, mutationTemplate.Sequence[i]));
            }
            //Create the mutations; use a list to store the newly mutated nucleotides

            List<PositionNucleotidePair> newSequence = new List<PositionNucleotidePair>(mutationTemplate.Sequence.Length);
            //int mutations = (int)Math.Round((((100.0 - (double)similarity) / 100.0) * (double)characters.Count), 0);
            for (int i = 0; i < mutationCount; i++)
            {
                int positionOfMutation = rand.Next(0, characters.Count);
                PositionNucleotidePair toMutate = characters[positionOfMutation];
                characters.RemoveAt(positionOfMutation);

                toMutate.nucleotide = GetTransversion(toMutate.nucleotide);
                newSequence.Add(toMutate);
            }


            //Write all other characters into the list containing the new sequence
            for (int i = 0; i < characters.Count; i++)
            {
                newSequence.Add(characters[i]);
            }



            //Sorter using anonymous delegate to sort 
            newSequence.Sort(delegate(PositionNucleotidePair p1, PositionNucleotidePair p2)
        {
            if (p1.position > p2.position) return 1;
            else if (p2.position > p1.position) return -1;
            else return 0;
        }
                );

            //Create a new ISequence
            ISequenceContainer retSeq = SequenceFactory.GetDefaultSequence();
            for (int i = 0; i < newSequence.Count; i++)
            {
                retSeq.Append(newSequence[i].nucleotide);
            }

            return new NucleotideSequence(mutationTemplate.Name + "_Mutated_Transversions;", retSeq);
        }


        #endregion
    }


    /// <summary>
    /// Utility class for trimming of NucleotideSequences,
    /// QualitySequences if available are also trimmed
    /// eg removal of polyA-tails etc
    /// </summary>
    public class TrimNucleotideSequence
    {


        public static NucleotideSequence Trim_5Prime_and_3PrimePolyA(NucleotideSequence toTrim, int toTrim5Prime, int maxDistanceFrom3PrimeEnd, int polyASeedSize)
        {
            ISequenceContainer seq = toTrim.Sequence;
            int toTrim3Prime = PolyA_ToTrim3Prime(toTrim.Sequence, maxDistanceFrom3PrimeEnd, polyASeedSize, polyASeedSize / 2);

            if (toTrim3Prime + toTrim5Prime >= toTrim.Length) return null;

            //123
            //ATCT
            //--
            int? startPos = null;
            if (toTrim.Start != null) startPos = toTrim.Start + toTrim5Prime;
            int? endPos = null;
            if (toTrim.End != null) endPos = toTrim.End - toTrim3Prime;
            NucleotideSequence toRet = new NucleotideSequence(toTrim.Name, toTrim.ParentName, toTrim.Comment, toTrim.Sequence.SubSequence(toTrim5Prime, seq.Length - toTrim5Prime - toTrim3Prime), startPos, endPos);

            //If the sequence contains a quality sequence trim this one also
            if (toTrim.QualitySequence != null)
            {
                QualitySequence oldQual = toTrim.QualitySequence;
                toRet.QualitySequence = new QualitySequence(oldQual.Name, oldQual.Sequence.SubSequence(toTrim5Prime, seq.Length - toTrim5Prime - toTrim3Prime));
            }

            return toRet;




        }


        public static List<NucleotideSequence> Trim_5Prime_and_3PrimePolyA(List<NucleotideSequence> toTrim, int toTrim5Prime, int maxDistanceFrom3PrimeEnd, int polyASeedSize)
        {
            List<NucleotideSequence> ret = new List<NucleotideSequence>(toTrim.Count);
            for (int i = 0; i < toTrim.Count; i++)
            {
                NucleotideSequence toAdd = Trim_5Prime_and_3PrimePolyA(toTrim[i], toTrim5Prime, maxDistanceFrom3PrimeEnd, polyASeedSize);
                if (toAdd != null) ret.Add(toAdd);
            }
            return ret;
        }


        private static int PolyA_ToTrim3Prime(ISequenceContainer seq, int maxDistanceFrom3PrimeEnd, int polyASeedSize, int maxDistanceBetweenTwoConsecutiveMismatches)
        {
            int? polyAStart = null;
            int lastMismatch = seq.Length + maxDistanceBetweenTwoConsecutiveMismatches + 1;
            bool seeded = false;
            int seedCount = 0;

            for (int i = seq.Length - 1; i >= 0; i--)
            {
                if (!seeded)
                {
                    if (seq[i] == 'A')//if the character is an A increase the count
                    {
                        seedCount++;
                        if (seedCount >= polyASeedSize) seeded = true; // if the 
                    }
                    else //if not set the A count to zero
                    {
                        seedCount = 0;
                        if (seq.Length - i > maxDistanceFrom3PrimeEnd)
                        {

                            goto end;
                        }

                    }


                }
                else
                {
                    if (seq[i] != 'A')
                    {
                        if (lastMismatch - i < maxDistanceBetweenTwoConsecutiveMismatches) //If the distance to the last mismatche is smaller than the distance specified the last mismatch is the end of the polyAtail
                        {
                            polyAStart = lastMismatch + 1;
                            goto end;
                        }     //01234
                        else////TAAAT
                        {
                            lastMismatch = i;
                        }
                    }
                }

            }

        end: ;

            if (polyAStart == null) return 0;
            //0123456
            //TCGTCGAAA //leng=9
            else return seq.Length - polyAStart.Value;

        }



        /// <summary>
        /// Trim a defined amount of base pairs from the 5' and from the 3' end of a nucleotide sequence
        /// </summary>
        /// <param name="toTrim">the sequence which should be trimmed</param>
        /// <param name="toTrim5Prime"></param>
        /// <param name="toTrim3Prime"></param>
        /// <returns></returns>
        public static NucleotideSequence Trim_5Prime_and_3Prime(NucleotideSequence toTrim, int toTrim5Prime, int toTrim3Prime)
        {
            ISequenceContainer seq = toTrim.Sequence;
            if (toTrim3Prime + toTrim5Prime >= toTrim.Length) return null;

            //123
            //ATCT
            //--
            int? startPos = null;
            if (toTrim.Start != null) startPos = toTrim.Start + toTrim5Prime;
            int? endPos = null;
            if (toTrim.End != null) endPos = toTrim.End - toTrim3Prime;
            NucleotideSequence toRet = new NucleotideSequence(toTrim.Name, toTrim.ParentName, toTrim.Comment, toTrim.Sequence.SubSequence(toTrim5Prime, seq.Length - toTrim5Prime - toTrim3Prime), startPos, endPos);
            //If a quality sequence has been assigned also trim the quality sequence
            if (toTrim.QualitySequence != null)
            {
                QualitySequence oldQual = toTrim.QualitySequence;
                toRet.QualitySequence = new QualitySequence(oldQual.Name, oldQual.Sequence.SubSequence(toTrim5Prime, seq.Length - toTrim5Prime - toTrim3Prime));
            }
            return toRet;
        }

        public static List<NucleotideSequence> Trim_5Prime_and_3Prime(List<NucleotideSequence> toTrim, int toTrim5Prime, int toTrim3Prime)
        {
            List<NucleotideSequence> ret = new List<NucleotideSequence>(toTrim.Count);
            for (int i = 0; i < toTrim.Count; i++)
            {
                ret.Add(Trim_5Prime_and_3Prime(toTrim[i], toTrim5Prime, toTrim3Prime));
            }
            return ret;
        }



    }





}