package casekit.nmr.dbservice;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Utils {

    public static int splitSDFile(final String pathToSDFile, final int maxMolPerFile) throws IOException, CDKException {
        final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToSDFile),
                SilentChemObjectBuilder.getInstance());
        final String fileEnding = pathToSDFile.split("\\.")[1];
        BufferedWriter bufferedWriter;
        SDFWriter sdfWriter;
        int counter = 0;
        int part = 1;
        bufferedWriter = new BufferedWriter(new FileWriter(pathToSDFile.split("\\.")[0]
                + "_"
                + part
                + "."
                + fileEnding));
        sdfWriter = new SDFWriter(bufferedWriter);
        while (iterator.hasNext()) {
            if (counter >= maxMolPerFile) {
                sdfWriter.close();

                part++;
                bufferedWriter = new BufferedWriter(new FileWriter(pathToSDFile.split("\\.")[0]
                        + "_"
                        + part
                        + "."
                        + fileEnding));
                sdfWriter = new SDFWriter(bufferedWriter);
                sdfWriter.write(iterator.next());
                counter = 1;
            } else {
                sdfWriter.write(iterator.next());
            }
            counter++;
        }
        sdfWriter.close();
        iterator.close();

        return part
                - 1;
    }
}
