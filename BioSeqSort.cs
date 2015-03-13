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

namespace Bio.Seq.Sort
{


    /// <summary>
    /// Sort a collection of positionables according the length of the feature
    /// -> Descending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_Length_Descending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {
            int lengthA = seqA.End.Value - seqA.Start.Value + 1;
            int lengthB = seqB.End.Value - seqB.Start.Value + 1;
            if (lengthA < lengthB) return 1;
            else if (lengthA > lengthB) return -1;
            else return 0;
        }
    }


    /// <summary>
    /// Sort a collection of positionables according the length of the feature
    /// -> Ascending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_Length_Ascending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {
            int lengthA = seqA.End.Value - seqA.Start.Value + 1;
            int lengthB = seqB.End.Value - seqB.Start.Value + 1;
            if (lengthA > lengthB) return 1;
            else if (lengthA < lengthB) return -1;
            else return 0;
        }
    }




    /// <summary>
    /// Sort a collection of items which implement the IPositionable interface
    /// 1..According to the name of the parent sequence
    /// 2..According to the start position of the positionable in the parent sequence
    /// -> Ascending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_ParentName_StartPos_Ascending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {

            int stringcompare = seqA.ParentName.CompareTo(seqB.ParentName);

            if (stringcompare != 0) return stringcompare;
            else
            {

                if (seqA.Start < seqB.Start) return -1;
                else if (seqA.Start > seqB.Start) return 1;
                else return 0;
            }
        }
    }

    /// <summary>
    /// Sort a collection of items which implement the IPositionable interface
    /// 1..According to the name of the parent sequence
    /// 2..According to the start position of the positionable in the parent sequence 
    /// -> Descending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_ParentName_StartPos_Descending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {

            int stringcompare = -seqA.ParentName.CompareTo(seqB.ParentName);

            if (stringcompare != 0) return stringcompare;
            else
            {

                if (seqA.Start < seqB.Start) return 1;
                else if (seqA.Start > seqB.Start) return -1;
                else return 0;
            }
        }
    }


    /*
    /// <summary>
    /// Sort a collection of items which implement the IPositionable interface
    /// 1..According to the name of the parent sequence
    /// 2..According to the start position of the positionable in the parent sequence - Ascending
    /// 3..According to the indel shift
    /// -> Ascending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_ParentName_StartPos_IndelShift_Ascending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {

            int stringcompare = seqA.ParentName.CompareTo(seqB.ParentName);

            if (stringcompare != 0) return stringcompare;
            else
            {

                if (seqA.Start < seqB.Start) return -1;
                else if (seqA.Start > seqB.Start) return 1;
                else
                {
                    if (seqA.Start_IndelShift < seqB.Start_IndelShift) return -1;
                    else if (seqA.Start_IndelShift > seqB.Start_IndelShift) return 1;
                    return 0;

                }
                    

            }
        }
    }

    /// <summary>
    /// Sort a collection of items which implement the IPositionable interface
    /// 1..According to the name of the parent sequence
    /// 2..According to the start position of the positionable in the parent sequence
    /// 3..According to the indelshift
    /// -> - Descending
    /// </summary>
    /// <typeparam name="TIPos"></typeparam>
    public class SortPositionables_ParentName_StartPos_IndelShift_Descending<TIPos> : IComparer<TIPos> where TIPos : IPositionable
    {
        int IComparer<TIPos>.Compare(TIPos seqA, TIPos seqB)
        {

            int stringcompare = -seqA.ParentName.CompareTo(seqB.ParentName);

            if (stringcompare != 0) return stringcompare;
            else
            {

                if (seqA.Start < seqB.Start) return 1;
                else if (seqA.Start > seqB.Start) return -1;
                else{
                    if (seqA.Start_IndelShift<seqB.Start_IndelShift) return 1;
                    else if(seqA.Start_IndelShift>seqB.Start_IndelShift) return -1;
                    else return 0;


                }
            }
        }
    }
     * 
     */





}