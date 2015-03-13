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
using Bio.Seq;
using Bio.Alignment;

namespace Bio.SNP.Sort
{


    /// <summary>
    /// Sort the sibbling SNPs according to the frequency of all
    /// </summary>
    public class SortSibblingSNPs_SNPFrequencyAll : IComparer<SibblingSNPs>
    {
        public SortSibblingSNPs_SNPFrequencyAll()
        {
        }

        int IComparer<SibblingSNPs>.Compare(SibblingSNPs ss1, SibblingSNPs ss2)
        {
            if (ss1.Count_SNPs_Valid < ss2.Count_SNPs_Valid)
            {
                return 1;
            }
            else if (ss1.Count_SNPs_Valid > ss2.Count_SNPs_Valid)
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
    }



    /// <summary>
    /// Sort the transcripts according to the frequency in the subset
    /// </summary>
    public class SortSNPSites_GeneID_Frequency : IComparer<SNPSite>
    {
        public SortSNPSites_GeneID_Frequency()
        {
        }

        int IComparer<SNPSite>.Compare(SNPSite sl1, SNPSite sl2)
        {

            ///Proceed in a threefold manner
            ///1...sort according to the gene ID
            ///2...sort according to the frequency at this locus in the subset
            int comp = String.Compare(sl1.DatabaseParent, sl2.DatabaseParent);

            if (comp == 0)
            {

                if (sl1.Count_Valid < sl2.Count_Valid)
                {
                    return 1;
                }
                else if (sl1.Count_Valid > sl2.Count_Valid)
                {
                    return -1;
                }
                else
                {
                    return 0;
                }
            }
            else
            {

                return comp;

            }
        }
    }

    /// <summary>
    /// Sort the transcripts according to the PIC values
    /// </summary>
    public class SortSNPSites_GeneID_PIC : IComparer<SNPSite>
    {
        public SortSNPSites_GeneID_PIC()
        {
        }

        int IComparer<SNPSite>.Compare(SNPSite sl1, SNPSite sl2)
        {

            ///Proceed in a threefold manner
            ///1...sort according to the gene ID
            ///2...sort according to the frequency at this locus in the subset
            int comp = String.Compare(sl1.DatabaseParent, sl2.DatabaseParent);

            if (comp == 0)
            {

                if (sl1.PIC < sl2.PIC)
                {
                    return 1;
                }
                else if (sl1.PIC > sl2.PIC)
                {
                    return -1;
                }
                else
                {
                    return 0;
                }

            }
            else
            {

                return comp;

            }
        }
    }


    /// <summary>
    /// Sort the transcripts according to the PIC values
    /// </summary>
    public class SortSNPSites_PIC : IComparer<SNPSite>
    {
        public SortSNPSites_PIC()
        {
        }

        int IComparer<SNPSite>.Compare(SNPSite sl1, SNPSite sl2)
        {


            if (sl1.PIC < sl2.PIC)
            {
                return 1;
            }
            else if (sl1.PIC > sl2.PIC)
            {
                return -1;
            }
            else
            {
                return 0;
            }
        }
    }

    /// <summary>
    /// Sort the transcripts according to the frequency in the subset
    /// </summary>
    public class SortSNPSite_Frequency : IComparer<SNPSite>
    {
        public SortSNPSite_Frequency()
        {
        }

        int IComparer<SNPSite>.Compare(SNPSite s1, SNPSite s2)
        {



            if (s1.Count_Valid < s2.Count_Valid)
            {
                return 1;
            }
            else if (s1.Count_Valid > s2.Count_Valid)
            {
                return -1;
            }
            else
            {
                return 0;
            }


        }
    }


    /// <summary>
    /// Sort the transcripts according to the frequency in the subset
    /// </summary>
    public class SortSNPs_DatabaseName_DatabasePosition_Allele : IComparer<SNP>
    {
        public SortSNPs_DatabaseName_DatabasePosition_Allele()
        {
        }

        int IComparer<SNP>.Compare(SNP snp1, SNP snp2)
        {

            ///Proceed in a threefold manner
            ///1...sort according to the database parent
            ///2...sort according to the start position of the database
            ///3...sort according to the snp character = allele
            int comp = String.Compare(snp1.DatabaseParent, snp2.DatabaseParent);

            if (comp == 0)
            {

                if (snp1.SiteBasis_Database < snp2.SiteBasis_Database)
                {
                    return -1;
                }
                else if (snp1.SiteBasis_Database > snp2.SiteBasis_Database)
                {
                    return 1;
                }
                else
                {
                    if (snp1.SiteIndelShift_Database < snp2.SiteIndelShift_Database)
                    {
                        return -1;
                    }
                    else if (snp1.SiteIndelShift_Database > snp2.SiteIndelShift_Database)
                    {
                        return 1;
                    }
                    else
                    {

                        if ((byte)snp1.QueryChar > (byte)snp2.QueryChar) return 1;
                        else if ((byte)snp1.QueryChar < (byte)snp2.QueryChar) return -1;
                        else return 0;
                    }
                }

            }
            else
            {

                return comp;

            }
        }
    }

    /// <summary>
    /// Sort the transcripts according to the frequency in the subset
    /// </summary>
    public class SortSNPs_DatabaseName_DatabasePosition : IComparer<SNP>
    {
        public SortSNPs_DatabaseName_DatabasePosition()
        {
        }

        int IComparer<SNP>.Compare(SNP snp1, SNP snp2)
        {

            ///Proceed in a threefold manner
            ///1...sort according to the database parent
            ///2...sort according to the start position of the database
            ///3...sort according to the snp character = allele
            int comp = String.Compare(snp1.DatabaseParent, snp2.DatabaseParent);

            if (comp == 0)
            {

                if (snp1.SiteBasis_Database < snp2.SiteBasis_Database)
                {
                    return -1;
                }
                else if (snp1.SiteBasis_Database > snp2.SiteBasis_Database)
                {
                    return 1;
                }
                else
                {
                    if (snp1.SiteIndelShift_Database < snp2.SiteIndelShift_Database)
                    {
                        return -1;
                    }
                    else if (snp1.SiteIndelShift_Database > snp2.SiteIndelShift_Database)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                }

            }
            else
            {

                return comp;

            }
        }
    }


    /// <summary>
    /// Sort the transcripts according to the frequency in the subset
    /// </summary>
    public class SortSNPList_Count : IComparer<List<SNP>>
    {
        public SortSNPList_Count()
        {
        }

        int IComparer<List<SNP>>.Compare(List<SNP> li1, List<SNP> li2)
        {

            if (li1.Count < li2.Count)
            {
                return 1;
            }
            else if (li1.Count > li2.Count)
            {
                return -1;
            }
            else
            {

                return 0;
            }


        }
    }
}