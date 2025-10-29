/*
 * The MIT License
 *
 * Copyright 2019 Michael Wenk [https://github.com/michaelwenk]
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package casekit.nmr.dbservice;

import casekit.nmr.model.*;
import casekit.nmr.utils.Utils;

import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class NMRShiftDB {

    public static String getSolvent(final String solventPropertyString, final String spectrumIndexInRecord) {
        final String[] solventPropertyStringSplit = solventPropertyString.split(":");
        String solvent;
        for (int i = 0; i < solventPropertyStringSplit.length - 1; i++) {
            if (solventPropertyStringSplit[i].trim().endsWith(spectrumIndexInRecord)) {
                solvent = solventPropertyStringSplit[i + 1].trim();
                while (solvent.substring(solvent.length()
                        - 1)
                        .matches("\\d")) {
                    solvent = solvent.substring(0, solvent.length()
                            - 1);
                }
                solvent = solvent.trim();

                return solvent;
            }
        }

        return null;
    }

    public static List<String> getSpectraProperties1D(final IAtomContainer ac, final String nucleus) {
        final List<String> spectraProperties1D = new ArrayList<>();
        for (final Object obj : ac.getProperties()
                .keySet()) {
            if (obj instanceof String
                    && ((String) obj).startsWith("Spectrum "
                            + nucleus)) {
                spectraProperties1D.add((String) obj);
            }
        }

        return spectraProperties1D;
    }

    /**
     * Returns a {@link DataSet} class object
     * for each valid molecule record in the given NMRShiftDB file. Valid means
     * here that each molecule record has to contain the given spectrum
     * property string as well as the number of signals in that spectrum has to
     * be the same as atoms of that atom type in molecule.
     *
     * @param pathToNMRShiftDB path to NMRShiftDB file
     * @param nuclei           nuclei to get the spectra for
     *
     * @return
     *
     * @throws CDKException
     * @throws IOException
     * @see DataSet
     */
    public static List<DataSet> getDataSetsFromNMRShiftDB(final String pathToNMRShiftDB,
            final String[] nuclei) throws CDKException, IOException {
        final List<DataSet> dataSets = new ArrayList<>();
        final IteratingSDFReader iterator = new IteratingSDFReader(new FileReader(pathToNMRShiftDB),
                SilentChemObjectBuilder.getInstance());
        IAtomContainer structure;
        Spectrum spectrum;
        Assignment assignment;
        DataSet dataSet;
        List<String> spectraProperties1D;
        String[] split;
        String spectrumIndexInRecord;
        List<Integer> explicitHydrogenIndices;
        int[] temp;

        while (iterator.hasNext()) {
            structure = iterator.next();
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(structure);
            explicitHydrogenIndices = casekit.nmr.utils.Utils.getExplicitHydrogenIndices(structure);
            Collections.sort(explicitHydrogenIndices);
            dataSet = Utils.atomContainerToDataSet(structure);

            for (final String nucleus : nuclei) {
                spectraProperties1D = getSpectraProperties1D(structure, nucleus);
                for (final String spectrumProperty1D : spectraProperties1D) {
                    split = spectrumProperty1D.split("\\s");
                    spectrumIndexInRecord = split[split.length
                            - 1];

                    // skip molecules which do not contain any of requested spectrum information
                    spectrum = NMRShiftDBSpectrumToSpectrum(structure.getProperty(spectrumProperty1D), nucleus);
                    // if no spectrum could be built or the number of signals in spectrum is
                    // different than the atom number in molecule
                    if ((spectrum == null)
                            || casekit.nmr.utils.Utils.getDifferenceSpectrumSizeAndMolecularFormulaCount(spectrum,
                                    Utils.getMolecularFormulaFromString(
                                            dataSet.getMeta()
                                                    .get("mf")),
                                    0) != 0) {
                        continue;
                    }
                    if (structure.getProperty("Solvent") != null) {
                        spectrum.addMetaInfo("solvent",
                                getSolvent(structure.getProperty("Solvent"), spectrumIndexInRecord));
                    }
                    if (structure.getProperty("Field Strength [MHz]") != null) {
                        for (final String fieldStrength : structure.getProperty("Field Strength [MHz]")
                                .toString()
                                .split("\\s")) {
                            if (fieldStrength.startsWith(spectrumIndexInRecord
                                    + ":")) {
                                try {
                                    spectrum.addMetaInfo("spectrometerFrequency", fieldStrength.split(
                                            spectrumIndexInRecord
                                                    + ":")[1]);
                                } catch (final NumberFormatException e) {
                                    // e.printStackTrace();
                                }
                                break;
                            }
                        }
                    }

                    assignment = NMRShiftDBSpectrumToAssignment(structure.getProperty(spectrumProperty1D), nucleus);
                    if (assignment != null
                            && !explicitHydrogenIndices.isEmpty()) {
                        int hCount;
                        for (int i = 0; i < assignment.getSize(); i++) {
                            for (int k = 0; k < assignment.getAssignment(0, i).length; k++) {
                                hCount = 0;
                                for (int j = 0; j < explicitHydrogenIndices.size(); j++) {
                                    if (explicitHydrogenIndices.get(j) >= assignment.getAssignment(0, i, k)) {
                                        break;
                                    }
                                    hCount++;
                                }
                                temp = assignment.getAssignment(0, i);
                                temp[k] = assignment.getAssignment(0, i, k)
                                        - hCount;
                                assignment.setAssignment(0, i, temp);
                            }
                        }
                    }
                    dataSet.setSpectrum(new SpectrumCompact(spectrum));
                    dataSet.setAssignment(assignment);

                    dataSets.add(dataSet.buildClone());
                }
            }
        }

        iterator.close();

        return dataSets;
    }

    // /**
    // * Returns a hashmap containing combined keys (by "_") of solvents
    // * and lists of calculated deviations between all given spectra for a
    // * nucleus in molecule record as values. <br>
    // * Here, only molecule records in NMRShiftDB file are considered which have
    // * at least two different spectra for same nucleus. <br>
    // * Example: "Spectrum 13C 0", "Spectrum 13C 1" will be used for given
    // * nucleus 13C.
    // *
    // * @param pathToNMRShiftDB
    // * @param nucleus
    // *
    // * @return
    // *
    // * @throws FileNotFoundException
    // * @throws CDKException
    // */
    // public static HashMap<String, ArrayList<Double>> getSolventDeviations(final
    // String pathToNMRShiftDB, final String nucleus) throws FileNotFoundException,
    // CDKException {
    // int signalCount;
    // Spectrum spectrum;
    // Assignment assignment;
    // final ArrayList<ArrayList<Object[]>> spectraSets =
    // getSpectraFromNMRShiftDB(pathToNMRShiftDB, nucleus);
    // HashMap<Integer, ArrayList<Double>> shiftsPerAtom;
    // HashMap<Integer, ArrayList<String>> solventsPerAtom;
    // ArrayList<String> solvents;
    // String[] solventsToSort;
    //
    // final HashMap<String, ArrayList<Double>> deviations = new HashMap<>();
    // String combiKey;
    //
    // for (final ArrayList<Object[]> spectraSetInRecord : spectraSets) {
    // shiftsPerAtom = new HashMap<>();
    // solventsPerAtom = new HashMap<>();
    // signalCount = -1;
    // for (final Object[] spectrumAndAssignment : spectraSetInRecord) {
    // spectrum = (Spectrum) spectrumAndAssignment[0];
    // assignment = (Assignment) spectrumAndAssignment[1];
    // if (signalCount == -1) {
    // signalCount = spectrum.getSignalCount();
    // } else if (signalCount != spectrum.getSignalCount()) {
    // continue;
    // }
    // for (final int atomIndex : assignment.getAssignments(0)) {
    // if (!shiftsPerAtom.containsKey(atomIndex)) {
    // shiftsPerAtom.put(atomIndex, new ArrayList<>());
    // solventsPerAtom.put(atomIndex, new ArrayList<>());
    // }
    // shiftsPerAtom.get(atomIndex).add(spectrum.getSignal(assignment.getIndex(0,
    // atomIndex)).getShift(0));
    // solventsPerAtom.get(atomIndex).add(spectrum.getSolvent());
    // }
    // }
    // if (shiftsPerAtom.isEmpty() ||
    // (shiftsPerAtom.get(Collections.min(shiftsPerAtom.keySet())).size() < 2)) {
    // continue;
    // }
    // solvents = new
    // ArrayList<>(solventsPerAtom.get(Collections.min(solventsPerAtom.keySet())));
    // // if(Collections.frequency(solvents, "Unreported") +
    // Collections.frequency(solvents, "Unknown") > solvents.size() - 2){
    // // continue;
    // // }
    //
    // for (final int atomIndex : shiftsPerAtom.keySet()) {
    // for (int s1 = 0; s1 < solvents.size(); s1++) {
    // // if(solvents.get(s1).equals("Unreported") ||
    // solvents.get(s1).equals("Unknown")){
    // // continue;
    // // }
    // for (int s2 = s1 + 1; s2 < solvents.size(); s2++) {
    // // if (solvents.get(s2).equals("Unreported") ||
    // solvents.get(s2).equals("Unknown")) {
    // // continue;
    // // }
    // solventsToSort = new String[2];
    // solventsToSort[0] = solvents.get(s1);
    // solventsToSort[1] = solvents.get(s2);
    // Arrays.sort(solventsToSort);
    // combiKey = solventsToSort[0] + "_" + solventsToSort[1];
    // if (!deviations.containsKey(combiKey)) {
    // deviations.put(combiKey, new ArrayList<>());
    // }
    // deviations.get(combiKey).add(Math.abs(shiftsPerAtom.get(atomIndex).get(s1) -
    // shiftsPerAtom.get(atomIndex).get(s2)));
    // }
    // }
    // }
    // }
    //
    // return deviations;
    // }
    //
    // /**
    // * @param pathToDB
    // *
    // * @return
    // *
    // * @throws FileNotFoundException
    // * @deprecated
    // */
    // public static Set<String> getAtomTypesInDB(final String pathToDB) throws
    // FileNotFoundException {
    // final HashSet<String> atomTypes = new HashSet<>();
    // final IteratingSDFReader iterator = new IteratingSDFReader(new
    // FileReader(pathToDB), SilentChemObjectBuilder.getInstance());
    // while (iterator.hasNext()) {
    // atomTypes.addAll(HOSECodeUtilities.getAtomTypesInAtomContainer(iterator.next()));
    // }
    //
    // return atomTypes;
    // }

    /**
     * Creates a two dimensional array of a given NMRShiftDB casekit.nmr entry
     * with all signal shift values, intensities, multiplicities and atom indices.
     *
     * @param NMRShiftDBSpectrum
     *
     * @return two dimensional array:
     *         1. dimension: signal index (row);
     *         2. dimension: signal shift value (column 1), signal intensity (column
     *         2),
     *         signal multiplicity (column 3), atom index in structure (column 4)
     */
    public static String[][] parseNMRShiftDBSpectrum(final String NMRShiftDBSpectrum) {
        if (NMRShiftDBSpectrum.trim()
                .isEmpty()) {
            return new String[][] {};
        }
        String[] signalSplit;
        final String[] shiftsSplit = NMRShiftDBSpectrum.split("\\|");
        final String[][] values = new String[shiftsSplit.length][4];
        for (int i = 0; i < shiftsSplit.length; i++) {
            signalSplit = shiftsSplit[i].split(";");
            values[i][0] = signalSplit[0]; // shift value
            values[i][1] = signalSplit[1].toLowerCase()
                    .split("[a-z]")[0]; // intensity
            values[i][2] = signalSplit[1].split("\\d+\\.\\d+").length > 0
                    ? signalSplit[1].split("\\d+\\.\\d+")[1].toLowerCase()
                    : ""; // multiplicity
            values[i][3] = signalSplit[2]; // atom index
        }

        return values;
    }

    public static Spectrum NMRShiftDBSpectrumToSpectrum(final String NMRShiftDBSpectrum, final String nucleus) {
        if ((NMRShiftDBSpectrum == null)
                || NMRShiftDBSpectrum.trim()
                        .isEmpty()) {
            return null;
        }
        final String[][] spectrumStringArray = parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = new Spectrum();
        spectrum.setNuclei(new String[] { nucleus });
        spectrum.setSignals(new ArrayList<>());
        String multiplicity;
        Double shift, intensity;
        try {
            for (int i = 0; i < spectrumStringArray.length; i++) {
                shift = Double.parseDouble(spectrumStringArray[i][0]);
                intensity = Double.parseDouble(spectrumStringArray[i][1]);
                multiplicity = spectrumStringArray[i][2].trim()
                        .isEmpty()
                                ? null
                                : spectrumStringArray[i][2].trim()
                                        .toLowerCase();
                spectrum.addSignal(
                        new Signal(new String[] { nucleus }, new Double[] { shift }, multiplicity, "signal", intensity,
                                1, 0,
                                null, null));
            }
        } catch (final Exception e) {
            return null;
        }

        return spectrum;
    }

    public static Assignment NMRShiftDBSpectrumToAssignment(final String NMRShiftDBSpectrum, final String nucleus) {
        if ((NMRShiftDBSpectrum == null)
                || NMRShiftDBSpectrum.trim()
                        .isEmpty()) {
            return null;
        }
        final String[][] NMRShiftDBSpectrumStringArray = parseNMRShiftDBSpectrum(NMRShiftDBSpectrum);
        final Spectrum spectrum = NMRShiftDBSpectrumToSpectrum(NMRShiftDBSpectrum, nucleus);
        final Assignment assignment = new Assignment();
        assignment.setNuclei(spectrum.getNuclei());
        assignment.initAssignments(spectrum.getSignalCount());
        int signalIndex;
        String multiplicity;
        List<Integer> closestSignalList;
        for (int i = 0; i < NMRShiftDBSpectrumStringArray.length; i++) {
            // just to be sure that we take the right signal if equivalences are present
            closestSignalList = spectrum.pickByClosestShift(Double.parseDouble(NMRShiftDBSpectrumStringArray[i][0]), 0,
                    0.0);
            multiplicity = NMRShiftDBSpectrumStringArray[i][2].trim()
                    .isEmpty()
                            ? null
                            : NMRShiftDBSpectrumStringArray[i][2].trim()
                                    .toLowerCase();
            closestSignalList.retainAll(spectrum.pickByMultiplicity(multiplicity));
            signalIndex = closestSignalList.get(0);

            assignment.addAssignmentEquivalence(0, signalIndex, Integer.parseInt(NMRShiftDBSpectrumStringArray[i][3]));
        }

        return assignment;
    }
}
