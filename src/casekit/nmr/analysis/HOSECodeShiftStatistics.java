package casekit.nmr.analysis;

import casekit.nmr.dbservice.COCONUT;
import casekit.nmr.dbservice.NMRShiftDB;
import casekit.nmr.fragments.model.ConnectionTree;
import casekit.nmr.hose.HOSECodeBuilder;
import casekit.nmr.model.DataSet;
import casekit.nmr.model.Signal;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Statistics;
import casekit.nmr.utils.Utils;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.reflect.TypeToken;
import org.bson.Document;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.nmrshiftdb.util.ExtendedHOSECodeGenerator;

import java.io.*;
import java.util.*;

public class HOSECodeShiftStatistics {

    private final static Gson GSON = new GsonBuilder().setLenient()
            .create();
    private final static ExtendedHOSECodeGenerator extendedHOSECodeGenerator = new ExtendedHOSECodeGenerator();

    public static Map<String, Map<String, List<Double>>> collectHOSECodeShifts(final List<DataSet> dataSetList,
            final Integer maxSphere,
            final boolean use3D,
            final boolean withExplicitH) {
        return collectHOSECodeShifts(dataSetList, maxSphere, use3D, withExplicitH, new HashMap<>());
    }

    /**
     * This method expects datasets containing structures without explicit
     * hydrogens.
     *
     * @param dataSetList
     * @param maxSphere
     * @param hoseCodeShifts
     *
     * @return
     */
    public static Map<String, Map<String, List<Double>>> collectHOSECodeShifts(final List<DataSet> dataSetList,
            final Integer maxSphere,
            final boolean use3D,
            final boolean withExplicitH,
            final Map<String, Map<String, List<Double>>> hoseCodeShifts) {
        for (final DataSet dataSet : dataSetList) {
            insert(dataSet, maxSphere, use3D, withExplicitH, hoseCodeShifts);
        }

        return hoseCodeShifts;
    }

    public static boolean insert(final DataSet dataSet, final Integer maxSphere, final boolean use3D,
            final boolean withExplicitH,
            final Map<String, Map<String, List<Double>>> hoseCodeShifts) {
        final IAtomContainer structure;
        Signal signal;
        String hoseCode;
        final String atomTypeSpectrum;
        String solvent;
        final Map<Integer, Integer> atomIndexMap; // from explicit H to heavy atom
        ConnectionTree connectionTree;
        int maxSphereTemp;
        List<Integer> signalIndices;
        structure = dataSet.getStructure()
                .toAtomContainer();
        final Spectrum spectrum = dataSet.getSpectrum()
                .toSpectrum();
        if (Utils.containsExplicitHydrogens(structure)) {
            System.out.println("!!!Dataset skipped must not contain (previously set) explicit hydrogens!!!");
            return false;
        }
        // create atom index map to know which indices the explicit hydrogens will have
        atomIndexMap = new HashMap<>();
        if (use3D
                || withExplicitH) {
            try {
                int nextAtomIndexExplicitH = structure.getAtomCount();
                for (int i = 0; i < structure.getAtomCount(); i++) {
                    if (structure.getAtom(i)
                            .getImplicitHydrogenCount() != null) {
                        for (int j = 0; j < structure.getAtom(i)
                                .getImplicitHydrogenCount(); j++) {
                            atomIndexMap.put(nextAtomIndexExplicitH, i);
                            nextAtomIndexExplicitH++;
                        }
                    }
                }

                if (use3D) {
                    try {
                        Utils.placeExplicitHydrogens(structure);
                    } catch (final CDKException | IOException | ClassNotFoundException e) {
                        e.printStackTrace();
                    }
                } else {
                    Utils.convertImplicitToExplicitHydrogens(structure);
                }
                Utils.setAromaticityAndKekulize(structure);
            } catch (final CDKException e) {
                e.printStackTrace();
                return false;
            }
        }
        solvent = dataSet.getSpectrum()
                .getMeta() == null
                        ? null
                        : dataSet.getSpectrum()
                                .getMeta()
                                .get("solvent");
        if (solvent == null
                || solvent.equals("")) {
            solvent = "Unknown";
        }
        atomTypeSpectrum = Utils.getAtomTypeFromNucleus(dataSet.getSpectrum()
                .getNuclei()[0]);
        for (int i = 0; i < structure.getAtomCount(); i++) {
            signalIndices = null;
            if (structure.getAtom(i)
                    .getSymbol()
                    .equals(atomTypeSpectrum)) {
                if (atomTypeSpectrum.equals("H")) {
                    // could be multiple signals
                    signalIndices = dataSet.getAssignment()
                            .getIndices(0, atomIndexMap.get(i));
                } else {
                    // should be one only
                    signalIndices = dataSet.getAssignment()
                            .getIndices(0, i);
                }
            }
            if (signalIndices != null) {
                for (final Integer signalIndex : signalIndices) {
                    signal = spectrum.getSignal(signalIndex);
                    try {
                        if (maxSphere == null) {
                            connectionTree = HOSECodeBuilder.buildConnectionTree(structure, i, null);
                            maxSphereTemp = connectionTree.getMaxSphere(true);
                        } else {
                            maxSphereTemp = maxSphere;
                        }
                        for (int sphere = 1; sphere <= maxSphereTemp; sphere++) {
                            if (use3D) {
                                try {
                                    hoseCode = extendedHOSECodeGenerator.getHOSECode(structure, structure.getAtom(i),
                                            sphere);
                                } catch (final Exception e) {
                                    // e.printStackTrace();
                                    continue;
                                }
                            } else {
                                hoseCode = HOSECodeBuilder.buildHOSECode(structure, i, sphere, false);
                            }
                            hoseCodeShifts.putIfAbsent(hoseCode, new HashMap<>());
                            hoseCodeShifts.get(hoseCode)
                                    .putIfAbsent(solvent, new ArrayList<>());
                            hoseCodeShifts.get(hoseCode)
                                    .get(solvent)
                                    .add(signal.getShift(0));
                        }
                    } catch (final CDKException e) {
                        e.printStackTrace();
                    }
                }
            }
        }

        return true;
    }

    public static Map<String, Map<String, Double[]>> buildHOSECodeShiftStatistics(
            final Map<String, Map<String, List<Double>>> hoseCodeShifts) {

        final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        List<Double> values;
        for (final Map.Entry<String, Map<String, List<Double>>> hoseCodes : hoseCodeShifts.entrySet()) {
            hoseCodeShiftStatistics.put(hoseCodes.getKey(), new HashMap<>());
            for (final Map.Entry<String, List<Double>> solvents : hoseCodes.getValue()
                    .entrySet()) {
                values = new ArrayList<>(solvents.getValue());
                values = Statistics.removeOutliers(values, 1.5);
                hoseCodeShiftStatistics.get(hoseCodes.getKey())
                        .put(solvents.getKey(),
                                new Double[] { (double) values.size(), Collections.min(values),
                                        Statistics.getMean(values), Statistics.getMedian(values),
                                        Collections.max(values) });
            }
        }

        return hoseCodeShiftStatistics;
    }

    public static Map<String, Map<String, Double[]>> buildHOSECodeShiftStatistics(final String[] pathsToNMRShiftDBs,
            final String[] pathsToCOCONUTs,
            final String[] nuclei,
            final Integer maxSphere,
            final boolean use3D,
            final boolean withExplicitH) {
        try {
            final Map<String, Map<String, List<Double>>> hoseCodeShifts = new HashMap<>();
            for (final String pathsToNMRShiftDB : pathsToNMRShiftDBs) {
                HOSECodeShiftStatistics.collectHOSECodeShifts(
                        NMRShiftDB.getDataSetsFromNMRShiftDB(pathsToNMRShiftDB, nuclei), maxSphere, use3D,
                        withExplicitH, hoseCodeShifts);
            }
            for (final String pathsToCOCONUT : pathsToCOCONUTs) {
                HOSECodeShiftStatistics.collectHOSECodeShifts(
                        COCONUT.getDataSetsWithShiftPredictionFromCOCONUT(pathsToCOCONUT, nuclei), maxSphere, use3D,
                        withExplicitH, hoseCodeShifts);
            }
            return HOSECodeShiftStatistics.buildHOSECodeShiftStatistics(hoseCodeShifts);
        } catch (final FileNotFoundException | CDKException e) {
            e.printStackTrace();
        }

        return new HashMap<>();
    }

    public static Map<String, Map<String, Double[]>> buildHOSECodeShiftStatistics(final List<DataSet> dataSetList,
            final Integer maxSphere,
            final boolean use3D,
            final boolean withExplicitH) {
        return HOSECodeShiftStatistics.buildHOSECodeShiftStatistics(
                collectHOSECodeShifts(dataSetList, maxSphere, use3D, withExplicitH));
    }

    public static boolean writeHOSECodeShiftStatistics(final Map<String, Map<String, Double[]>> hoseCodeShifts,
            final String pathToJsonFile) {
        try {
            final BufferedWriter bw = new BufferedWriter(new FileWriter(pathToJsonFile));
            bw.append("{");
            bw.newLine();
            bw.flush();

            Document subDocument;
            String json;
            long counter = 0;
            for (final Map.Entry<String, Map<String, Double[]>> entry : hoseCodeShifts.entrySet()) {
                subDocument = new Document();
                subDocument.append("HOSECode", entry.getKey());
                subDocument.append("values", GSON.toJson(entry.getValue()));
                json = new Document(String.valueOf(counter), subDocument).toJson();
                bw.append(json, 1, json.length()
                        - 1);
                if (counter < hoseCodeShifts.size()
                        - 1) {
                    bw.append(",");
                }
                bw.newLine();
                bw.flush();

                counter++;
            }

            bw.append("}");
            bw.flush();
            bw.close();

            return true;
        } catch (final IOException e) {
            e.printStackTrace();
        }

        return false;
    }

    public static Map<String, Map<String, Double[]>> readHOSECodeShiftStatistics(
            final String pathToJsonFile) throws FileNotFoundException {
        final BufferedReader br = new BufferedReader(new FileReader(pathToJsonFile));
        final Map<String, Map<String, Double[]>> hoseCodeShiftStatistics = new HashMap<>();
        // add all task to do
        br.lines()
                .forEach(line -> {
                    if ((line.trim()
                            .length() > 1)
                            || (!line.trim()
                                    .startsWith("{")
                                    && !line.trim()
                                            .endsWith("}"))) {
                        final StringBuilder hoseCodeShiftsStatisticInJSON = new StringBuilder();
                        if (line.endsWith(",")) {
                            hoseCodeShiftsStatisticInJSON.append(line, 0, line.length()
                                    - 1);
                        } else {
                            hoseCodeShiftsStatisticInJSON.append(line);
                        }
                        final JsonObject jsonObject = JsonParser.parseString(hoseCodeShiftsStatisticInJSON.substring(
                                hoseCodeShiftsStatisticInJSON.toString()
                                        .indexOf("{")))
                                .getAsJsonObject();
                        hoseCodeShiftStatistics.put(jsonObject.get("HOSECode")
                                .getAsString(),
                                GSON.fromJson(jsonObject.get("values")
                                        .getAsString(),
                                        new TypeToken<Map<String, Double[]>>() {
                                        }.getType()));
                    }
                });

        return hoseCodeShiftStatistics;
    }
}
