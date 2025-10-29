package casekit.nmr.utils;

import casekit.nmr.model.DataSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.List;

public class Parser {

    public static List<DataSet> parseSDFile(final String pathToFile)
            throws CDKException, IOException {
        return parseSDFile(new FileReader(pathToFile));
    }

    public static List<DataSet> parseSDFileContent(final String fileContent) throws CDKException, IOException {
        final InputStream inputStream = new ByteArrayInputStream(fileContent.getBytes(StandardCharsets.UTF_8));
        return parseSDFile(new InputStreamReader(inputStream));
    }

    public static List<DataSet> parseSDFile(final Reader fileReader) throws CDKException, IOException {
        final List<DataSet> dataSetList = new ArrayList<>();

        final IteratingSDFReader iterator = new IteratingSDFReader(fileReader, SilentChemObjectBuilder.getInstance());
        while (iterator.hasNext()) {
            dataSetList.add(Utils.atomContainerToDataSet(iterator.next()));
        }
        iterator.close();

        return dataSetList;
    }

    public static List<DataSet> parseSmilesFile(final String pathToFile) throws FileNotFoundException {
        return parseSmilesFile(new FileReader(pathToFile));
    }

    public static List<DataSet> parseSmilesFileContent(final String fileContent) {
        final InputStream inputStream = new ByteArrayInputStream(fileContent.getBytes(StandardCharsets.UTF_8));
        return parseSmilesFile(new InputStreamReader(inputStream));
    }

    public static List<DataSet> parseSmilesFile(final Reader fileReader) {
        final List<DataSet> dataSetList = new ArrayList<>();
        final BufferedReader bufferedReader = new BufferedReader(fileReader);
        final SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        bufferedReader.lines()
                .forEach(smiles -> {
                    try {
                        dataSetList.add(Utils.atomContainerToDataSet(smilesParser.parseSmiles(smiles)));
                    } catch (final CDKException e) {
                        e.printStackTrace();
                    }
                });

        return dataSetList;
    }

    public static List<String> smilesFileToList(final String pathToFile) throws FileNotFoundException {
        return smilesFileToList(new FileReader(pathToFile));
    }

    public static List<String> smilesFileContentToList(final String fileContent) {
        final InputStream inputStream = new ByteArrayInputStream(fileContent.getBytes(StandardCharsets.UTF_8));
        return smilesFileToList(new InputStreamReader(inputStream));
    }

    public static List<String> smilesFileToList(final Reader fileReader) {
        final List<String> smilesList = new ArrayList<>();
        final BufferedReader bufferedReader = new BufferedReader(fileReader);
        bufferedReader.lines()
                .forEach(smilesList::add);

        return smilesList;
    }

}
