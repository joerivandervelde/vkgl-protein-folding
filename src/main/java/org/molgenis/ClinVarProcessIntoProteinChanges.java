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

        File clinvarBgzipAnnot = new File("/Users/joeri/ClinVar/clinvar_20230722_annot.vcf.gz");
        File outputFile = new File("/Users/joeri/ClinVar/clinvar_20230722_protForFolding.tsv");

        Map<String,String> protChangesToClsf = new TreeMap<>();
        Map<String,String> protChangesToChromPosRefAlt = new TreeMap<>();

//        InputStream fileStream = new FileInputStream(clinvarBgzipAnnot);
//        InputStream gzipStream = new GZIPInputStream(fileStream);
//        Reader decoder = new InputStreamReader(gzipStream, StandardCharsets.UTF_8);
//        BufferedReader buffered = new BufferedReader(decoder);

        BufferedReader in = new BufferedReader(new InputStreamReader(new BlockCompressedInputStream(clinvarBgzipAnnot)));


        String line;
        while ((line = in.readLine()) != null) {
        {
            if(line.startsWith("#"))
            {
                continue;
            }
            //line = s.nextLine();
            String[] lineSplit = line.split("\t", -1);
            if(lineSplit.length != 8)
            {
                throw new Exception("line split expecting 8 elements for " + line);
            }
            System.out.println(line);

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

}
