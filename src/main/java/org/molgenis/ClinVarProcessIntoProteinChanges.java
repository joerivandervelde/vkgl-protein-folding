package org.molgenis;

import htsjdk.samtools.util.BlockCompressedInputStream;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.Map;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;


public class ClinVarProcessIntoProteinChanges {

    public static void main(String args[]) throws Exception {
        System.out.println("Starting...");

        File clinvarBgzipAnnot = new File("/Users/joeri/ClinVar/clinvar_b38_20230722_vep.vcf.gz");
        File outputFile = new File("/Users/joeri/ClinVar/clinvar_20230722_protForFolding.tsv");

        Map<String,String> protChangesToClsf = new TreeMap<>();
        Map<String,String> protChangesToChromPosRefAlt = new TreeMap<>();

//        InputStream fileStream = new FileInputStream(clinvarBgzipAnnot);
//        InputStream gzipStream = new GZIPInputStream(fileStream);
//        Reader decoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
//        BufferedReader buffered = new BufferedReader(decoder);

        BufferedReader in = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(clinvarBgzipAnnot)));


        String line;
        loopOverLines:
        while ((line = in.readLine()) != null) {
            if(line.startsWith("#") ||
                    !line.contains("missense_variant") ||
                    line.contains("CLNSIG=Conflicting_interpretations_of_pathogenicity") ||
                    line.contains("CLNSIG=no_assertion_provided") ||
                    line.contains("CLNSIG=drug_response") ||
                    line.contains("CLNSIG=no_assertion_criteria_provided") ||
                    line.contains("CLNREVSTAT=no_interpretation_for_the_single_variant") ||
                    line.contains("CLNSIG=association") ||
                    line.contains("CLNSIG=Affects") ||
                    line.contains("CLNSIG=other") ||
                    line.contains("CLNSIG=not_provided") ||
                    line.contains("CLNSIG=risk_factor") ||
                    line.contains("CLNSIG=risk_allele") ||
                    line.contains("CLNSIG=Likely_risk_allele") ||
                    line.contains("CLNSIG=Uncertain_risk_allele") ||
                    line.contains("CLNSIG=confers_sensitivity") ||
                    line.contains("CLNSIG=protective")
            )
            {
                continue;
            }

            String[] lineSplit = line.split("\t", -1);
            if(lineSplit.length != 8)
            {
                throw new Exception("line split expecting 8 elements for " + line);
            }

            // length is between 12 and 21 or so, can't check on it
            String[] infoSplit = lineSplit[7].split(";", -1);

            for(String info : infoSplit)
            {
                if(info.startsWith("CSQ=")){

                    // before we deal with the consequences, get the classification
                    String clf = null;
                    for(String infoAgain : infoSplit) {
                        if (infoAgain.startsWith("CLNSIG=")) {
                            if(infoAgain.contains("Benign")){
                                clf = "LB";
                            }
                            else if(infoAgain.contains("Likely_benign")){
                                clf = "LB";
                            }
                            else if(infoAgain.contains("Uncertain_significance")){
                                clf = "VUS";
                            }
                            else if(infoAgain.contains("Likely_pathogenic")){
                                clf = "LP";
                            }
                            else if(infoAgain.contains("Pathogenic")){
                                clf = "LP";
                            }
                            else{
                                throw new Exception("Unsupported CLNSIG: " + infoAgain);
                            }
                        }
                    }
                    if(clf == null)
                    {
                        throw new Exception("No CLNSIG for line: " + line);
                    }

                    // now handle the consequence field
                    String[] csqSplit = info.split(",", -1);
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
                            throw new Exception("csq parts split expecting 25 elements, but found "+csqParts.length+" for " + line);
                        }

                        String proteinPos = csqParts[14];

                        // some variants are interpreted in a non-coding context in which case proteinPos is empty
                        if(proteinPos.isBlank()){
                            continue;
                        }

                        // synonymous variants have only 1 letter ('A') instead of 'A/B'
                        String[] proteinChangeSplit = csqParts[15].split("/");
                        if(proteinChangeSplit.length != 2)
                        {
                            continue;
                        }

                        // as per FoldX notation: WT residue, chain, residue number, mutant residue (e.g. "CA1490Y;")
                        String geneSymbol = csqParts[3];
                        String geneProt = geneSymbol + "\t" + proteinChangeSplit[0] + "A" + proteinPos + proteinChangeSplit[1];
                        //String clf = infoSplit[0].replace("VKGL=","");

                        // unless there is a different protein notation for a transcript, the notation is duplicate
                        // filter here by checking if gene-prot combination is already present
                        // at the same time check for any potential conflicting interpretations for this notation
                        if(protChangesToClsf.containsKey(geneProt))
                        {
                            String prevClf = protChangesToClsf.get(geneProt);
                            if(!prevClf.equals(clf))
                            {
                                System.out.println("Difference in classification for key '" + geneProt + "' at " + line);
                                if(prevClf.equals("LB"))
                                {
                                    System.out.println("Solution: overriding LB with " + clf);
                                    protChangesToClsf.put(geneProt, clf);
                                }else if(prevClf.equals("VUS") && clf.equals("LP"))
                                {
                                    System.out.println("Solution: overriding VUS with " + clf);
                                    protChangesToClsf.put(geneProt, clf);
                                }else if(prevClf.equals("LP"))
                                {
                                    System.out.println("Solution: keeping LP");
                                }
                                else{
                                    // old is VUS, new is LB
                                    System.out.println("Solution: keeping VUS");
                                }
                            }else{
                                // no difference in classification, so no need to add or do anything
                            }
                        }else{
                            // first time add, also store genome coordinates
                            protChangesToClsf.put(geneProt, clf);
                            protChangesToChromPosRefAlt.put(geneProt, (lineSplit[0] + "\t" + lineSplit[1] + "\t" + lineSplit[3] + "\t" + lineSplit[4]));
                        }
                    }

                    break;
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
            writer.write("GRCh38" + "\t" + protChangesToChromPosRefAlt.get(key) + "\t" + key + "\t" + protChangesToClsf.get(key) + System.lineSeparator());
        }
        writer.flush();
        writer.close();


    }

}
