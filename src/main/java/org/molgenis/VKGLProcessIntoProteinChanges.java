package org.molgenis;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

public class VKGLProcessIntoProteinChanges {

    public static void main(String args[]) throws Exception {
        System.out.println("Starting...");

        File vkglMissense = new File("/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_annot_missense.vcf");
        File outputFile = new File("/Users/joeri/VKGL/VKGL-prot/VKGL_apr2023_protForFolding.tsv");

        Map<String,String> protChangesToClsf = new TreeMap<>();
        Map<String,String> protChangesToChromPosRefAlt = new TreeMap<>();
        Scanner s = new Scanner(vkglMissense);
        String line;
        while(s.hasNextLine())
        {
            line = s.nextLine();
            String[] lineSplit = line.split("\t", -1);
            if(lineSplit.length != 8)
            {
                throw new Exception("line split expecting 8 elements for " + line);
            }

            String[] infoSplit = lineSplit[7].split(";", -1);
            if(infoSplit.length != 2)
            {
                throw new Exception("info split expecting 2 elements for " + line);
            }

            String[] csqSplit = infoSplit[1].split(",", -1);
            if(csqSplit.length == 0)
            {
                throw new Exception("csq split expecting 1 or more elements for " + line);
            }

            for(String csq : csqSplit)
            {
                String[] csqParts = csq.split("\\|", -1);
                // CSQ=T|missense_variant|MODERATE|CLIC2|1193|Transcript|NM_001289.6|protein_coding|6/6||||905|718|240|A/T|Gca/Aca|||-1||EntrezGene|||
                if(csqParts.length != 25)
                {
                    throw new Exception("csq parts split expecting 25 elements for " + line);
                }

                String proteinPos = csqParts[14];

                // some variants are interpreted in a non-coding context in which case proteinPos is empty
                if(proteinPos.isBlank()){
                    continue;
                }
                String[] proteinChangeSplit = csqParts[15].split("/");
                // synonymous variants have only 1 letter ('A') instead of 'A/B'
                if(proteinChangeSplit.length != 2)
                {
                    continue;
                }

                String geneSymbol = csqParts[3];
                // as per FoldX notation: WT residue, chain, residue number, mutant residue (e.g. "CA1490Y;")
                String geneProt = geneSymbol + "\t" + proteinChangeSplit[0] + "A" + proteinPos + proteinChangeSplit[1];
                String clf = infoSplit[0].replace("VKGL=","");

                // unless there is a different protein notation for a transcript, the notation is duplicate
                // filter here by checking if gene-prot combination is already present
                // at the same time check for any potential conflicting interpretations for this notation
                if(protChangesToClsf.containsKey(geneProt))
                {
                    String prevClf = protChangesToClsf.get(geneProt);
                    if(!prevClf.equals(clf))
                    {
                        System.out.println("Conflicting classification for key '" + geneProt + "' at " + line);
                        protChangesToClsf.put(geneProt, "CF");

//                        if(prevClf.equals("LB"))
//                        {
//                            System.out.println("Solution: overriding LB with " + clf);
//                            protChangesToClsf.put(geneProt, clf);
//                        }else if(prevClf.equals("VUS") && clf.equals("LP"))
//                        {
//                            System.out.println("Solution: overriding VUS with " + clf);
//                            protChangesToClsf.put(geneProt, clf);
//                        }else if(prevClf.equals("LP"))
//                        {
//                            System.out.println("Solution: keeping LP");
//                        }
//                        else{
//                            // old is VUS, new is LB
//                            System.out.println("Solution: keeping VUS");
//                        }
                    }else{
                        // no difference in classification, so no need to add or do anything
                    }
                }else{
                    // first time add, also store genome coordinates
                    protChangesToClsf.put(geneProt, clf);
                    protChangesToChromPosRefAlt.put(geneProt, (lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + "\t" + lineSplit[4]));
                }
            }
        }

        /*
        Write results, excluding blacklisted protein changes
         */
        PrintWriter writer = new PrintWriter(outputFile, "UTF-8");
        writer.write("Assembly" + "\t" + "Chrom" + "\t" + "Pos" + "\t" + "Ref" + "\t" + "Alt" + "\t" + "Gene" + "\t" + "ProtChange" + "\t" + "Classification" + System.lineSeparator());
        for(String key : protChangesToClsf.keySet())
        {
            writer.write("GRCh37" + "\t" + protChangesToChromPosRefAlt.get(key) + "\t" + key + "\t" + protChangesToClsf.get(key) + System.lineSeparator());
        }
        writer.flush();
        writer.close();

    }

}
