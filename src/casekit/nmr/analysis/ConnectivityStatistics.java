package casekit.nmr.analysis;

import casekit.nmr.model.DataSet;
import casekit.nmr.model.Spectrum;
import casekit.nmr.utils.Utils;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;

public class ConnectivityStatistics {

    /**
     * @param dataSet
     * @param atomType
     * @param occurrenceStatistics multiplicity -> hybridization -> shift (int) ->
     *                             elemental composition (mf) -> connected atom
     *                             symbol -> [#found, #notFound]
     */
    public static void buildOccurrenceStatistics(final DataSet dataSet, final String atomType,
            final Map<String, Map<String, Map<Integer, Map<String, Map<String, Integer[]>>>>> occurrenceStatistics) {
        final Spectrum spectrum = dataSet.getSpectrum()
                .toSpectrum();
        final IAtomContainer structure = dataSet.getStructure()
                .toAtomContainer();
        final List<String> elements = buildElements(structure);
        final String elementsString = buildElementsString(elements);

        int shift, atomIndex;
        IAtom atom;
        String multiplicity, hybridization;
        Set<String> found, notFound;
        for (int signalIndex = 0; signalIndex < spectrum.getSignalCount(); signalIndex++) {
            shift = spectrum.getShift(signalIndex, 0)
                    .intValue();
            for (int equivalenceIndex = 0; equivalenceIndex < dataSet.getAssignment()
                    .getAssignment(0, signalIndex).length; equivalenceIndex++) {
                atomIndex = dataSet.getAssignment()
                        .getAssignment(0, signalIndex, equivalenceIndex);
                atom = structure.getAtom(atomIndex);
                if (atom.getSymbol()
                        .equals(atomType)) {
                    multiplicity = Utils.getMultiplicityFromProtonsCount(atom.getImplicitHydrogenCount());
                    if (multiplicity == null || atom.getHybridization() == null) {
                        continue;
                    }
                    multiplicity = multiplicity.toLowerCase();
                    hybridization = atom.getHybridization()
                            .name();
                    occurrenceStatistics.putIfAbsent(multiplicity, new ConcurrentHashMap<>());
                    occurrenceStatistics.get(multiplicity)
                            .putIfAbsent(hybridization, new ConcurrentHashMap<>());
                    occurrenceStatistics.get(multiplicity)
                            .get(hybridization)
                            .putIfAbsent(shift, new ConcurrentHashMap<>());
                    occurrenceStatistics.get(multiplicity)
                            .get(hybridization)
                            .get(shift)
                            .putIfAbsent(elementsString, new ConcurrentHashMap<>());
                    // check for connected hetero atoms
                    found = new HashSet<>();
                    for (final IAtom connectedAtom : structure.getConnectedAtomsList(atom)) {
                        if (connectedAtom.getSymbol()
                                .equals("H")) {
                            continue;
                        }
                        found.add(connectedAtom.getSymbol());
                    }
                    for (final String connectedAtomType : found) {
                        occurrenceStatistics.get(multiplicity)
                                .get(hybridization)
                                .get(shift)
                                .get(elementsString)
                                .putIfAbsent(connectedAtomType, new Integer[] { 0, 0 });
                        occurrenceStatistics.get(multiplicity)
                                .get(hybridization)
                                .get(shift)
                                .get(elementsString)
                                .get(connectedAtomType)[0] = occurrenceStatistics.get(multiplicity)
                                        .get(hybridization)
                                        .get(shift)
                                        .get(elementsString)
                                        .get(connectedAtomType)[0]
                                        + 1;
                    }
                    notFound = new HashSet<>(elements);
                    notFound.removeAll(found);
                    for (final String notConnectedAtomType : notFound) {
                        occurrenceStatistics.get(multiplicity)
                                .get(hybridization)
                                .get(shift)
                                .get(elementsString)
                                .putIfAbsent(notConnectedAtomType, new Integer[] { 0, 0 });
                        occurrenceStatistics.get(multiplicity)
                                .get(hybridization)
                                .get(shift)
                                .get(elementsString)
                                .get(notConnectedAtomType)[1] = occurrenceStatistics.get(multiplicity)
                                        .get(hybridization)
                                        .get(shift)
                                        .get(elementsString)
                                        .get(notConnectedAtomType)[1]
                                        + 1;
                    }
                }
            }
        }
    }

    /**
     * @param structure            structure to build and add the statistics from
     * @param heavyAtomsStatistics elemental composition (mf) -> connected atom pair
     *                             -> #found
     */
    public static void buildHeavyAtomsStatistics(final IAtomContainer structure,
            final Map<String, Map<String, Integer>> heavyAtomsStatistics) {
        final List<String> elements = buildElements(structure);
        final String elementsString = buildElementsString(elements);
        heavyAtomsStatistics.putIfAbsent(elementsString, new HashMap<>());

        for (final String combination : buildCombinations(elements)) {
            heavyAtomsStatistics.get(elementsString)
                    .putIfAbsent(combination, 0);
        }

        IAtom atom;
        String atomPairKey;
        final Map<String, Integer> found = new HashMap<>();
        for (int i = 0; i < structure.getAtomCount(); i++) {
            atom = structure.getAtom(i);
            // check for connected hetero atoms
            for (final IAtom connectedAtom : structure.getConnectedAtomsList(atom)) {
                if (connectedAtom.getSymbol()
                        .equals("H")) {
                    continue;
                }
                atomPairKey = buildAtomPairString(atom.getSymbol(), connectedAtom.getSymbol());
                found.putIfAbsent(atomPairKey, 0);
                found.put(atomPairKey, found.get(atomPairKey)
                        + 1);
            }
        }
        for (final String connectedAtomPair : found.keySet()) {
            heavyAtomsStatistics.get(elementsString)
                    .put(connectedAtomPair, heavyAtomsStatistics.get(elementsString)
                            .get(connectedAtomPair)
                            + found.get(connectedAtomPair)
                                    / 2); // divided by two since we count for both bond partners in the previous loop
        }
    }

    public static List<String> buildElements(final IAtomContainer structure) {
        final String mf = Utils.molecularFormularToString(Utils.getMolecularFormulaFromAtomContainer(structure));
        return buildElements(mf);
    }

    public static List<String> buildElements(final String mf) {
        final List<String> elements = new ArrayList<>(Utils.getMolecularFormulaElementCounts(mf)
                .keySet());
        elements.remove("H");
        Collections.sort(elements);

        return elements;
    }

    public static String buildElementsString(final List<String> elements) {
        return String.join(",", elements);
    }

    public static String buildAtomPairString(final String atomType1, final String atomType2) {
        final List<String> atomPairList = new ArrayList<>();
        atomPairList.add(atomType1);
        atomPairList.add(atomType2);
        Collections.sort(atomPairList);

        return String.join("_", atomPairList);
    }

    public static Set<String> buildCombinations(final List<String> elements) {
        final Set<String> combinations = new HashSet<>();
        for (final String element1 : elements) {
            for (final String element2 : elements) {
                combinations.add(buildAtomPairString(element1, element2));
            }
        }

        return combinations;
    }

    // /**
    // * @param dataSetList
    // * @param nucleus
    // * @param connectivityStatistics multiplicity -> hybridization -> shift (int)
    // -> connected atom symbol -> connected atom hybridization -> connected atom
    // protons count -> occurrence
    // */
    // @Deprecated
    // public static void buildConnectivityStatistics(final List<DataSet>
    // dataSetList, final String nucleus,
    // final Map<String, Map<String, Map<Integer, Map<String, Map<String,
    // Map<Integer, Integer>>>>>> connectivityStatistics) {
    // final String atomType = Utils.getAtomTypeFromNucleus(nucleus);
    // for (final DataSet dataSet : dataSetList) {
    // if (!dataSet.getSpectrum()
    // .getNuclei()[0].equals(nucleus)) {
    // continue;
    // }
    // buildConnectivityStatistics(dataSet, atomType, connectivityStatistics);
    // }
    // }
    //
    // /**
    // * @param dataSet
    // * @param atomType
    // * @param connectivityStatistics multiplicity -> hybridization -> shift (int)
    // -> connected atom symbol -> connected atom hybridization -> connected atom
    // protons count -> occurrence
    // */
    // @Deprecated
    // public static void buildConnectivityStatistics(final DataSet dataSet, final
    // String atomType,
    // final Map<String, Map<String, Map<Integer, Map<String, Map<String,
    // Map<Integer, Integer>>>>>> connectivityStatistics) {
    // final IAtomContainer structure = dataSet.getStructure()
    // .toAtomContainer();
    // final Spectrum spectrum = dataSet.getSpectrum()
    // .toSpectrum();
    // int shift, atomIndex;
    // IAtom atom;
    // String multiplicity, hybridization, connectedAtomType,
    // connectedAtomHybridization;
    // for (int signalIndex = 0; signalIndex
    // < spectrum.getSignalCount(); signalIndex++) {
    // shift = spectrum.getShift(signalIndex, 0)
    // .intValue();
    // for (int equivalenceIndex = 0; equivalenceIndex
    // < dataSet.getAssignment()
    // .getAssignment(0, signalIndex).length; equivalenceIndex++) {
    // atomIndex = dataSet.getAssignment()
    // .getAssignment(0, signalIndex, equivalenceIndex);
    // atom = structure.getAtom(atomIndex);
    // if (atom.getSymbol()
    // .equals(atomType)) {
    // multiplicity =
    // Utils.getMultiplicityFromProtonsCount(atom.getImplicitHydrogenCount());
    // if (multiplicity
    // == null) {
    // continue;
    // }
    // multiplicity = multiplicity.toLowerCase();
    // hybridization = atom.getHybridization()
    // .name();
    // connectivityStatistics.putIfAbsent(multiplicity, new ConcurrentHashMap<>());
    // connectivityStatistics.get(multiplicity)
    // .putIfAbsent(hybridization, new ConcurrentHashMap<>());
    // // check for connected hetero atoms
    // for (final IAtom connectedAtom : structure.getConnectedAtomsList(atom)) {
    // if (connectedAtom.getSymbol()
    // .equals("H")) {
    // continue;
    // }
    // connectedAtomType = connectedAtom.getSymbol();
    // if (connectedAtom.getHybridization()
    // == null) {
    // continue;
    // }
    // connectedAtomHybridization = connectedAtom.getHybridization()
    // .name();
    // connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .putIfAbsent(shift, new ConcurrentHashMap<>());
    // connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .putIfAbsent(connectedAtomType, new ConcurrentHashMap<>());
    // connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .get(connectedAtomType)
    // .putIfAbsent(connectedAtomHybridization, new ConcurrentHashMap<>());
    // connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .get(connectedAtomType)
    // .get(connectedAtomHybridization)
    // .putIfAbsent(connectedAtom.getImplicitHydrogenCount(), 0);
    // connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .get(connectedAtomType)
    // .get(connectedAtomHybridization)
    // .put(connectedAtom.getImplicitHydrogenCount(), connectivityStatistics.get(
    // multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .get(connectedAtomType)
    // .get(connectedAtomHybridization)
    // .get(connectedAtom.getImplicitHydrogenCount())
    // + 1);
    // }
    // }
    // }
    // }
    // }
    //
    // /**
    // * @param connectivityStatistics multiplicity -> hybridization -> shift (int)
    // -> connected atom symbol -> connected atom hybridization -> connected atom
    // protons count -> occurrence
    // * @param multiplicity
    // * @param hybridization
    // * @param shift
    // * @param molecularFormulaElements
    // *
    // * @return
    // */
    // @Deprecated
    // public static Map<String, Map<String, Map<Integer, Integer>>>
    // extractConnectivities(
    // final Map<String, Map<String, Map<Integer, Map<String, Map<String,
    // Map<Integer, Integer>>>>>> connectivityStatistics,
    // final String multiplicity, final String hybridization, final int shift,
    // final Set<String> molecularFormulaElements) {
    // final Map<String, Map<String, Map<Integer, Integer>>> extractedConnectivities
    // = new HashMap<>();
    // if (connectivityStatistics.containsKey(multiplicity)
    // && connectivityStatistics.get(multiplicity)
    // .containsKey(hybridization)
    // && connectivityStatistics.get(multiplicity)
    // .get(hybridization)
    // .containsKey(shift)) {
    // for (final Map.Entry<String, Map<String, Map<Integer, Integer>>> entry :
    // connectivityStatistics.get(
    // multiplicity)
    // .get(hybridization)
    // .get(shift)
    // .entrySet()) {
    // if (molecularFormulaElements.contains(entry.getKey())) {
    // extractedConnectivities.put(entry.getKey(), entry.getValue());
    // }
    // }
    // }
    //
    // return extractedConnectivities;
    // }
    //
    // @Deprecated
    // public static Map<String, Map<Integer, Map<Integer, Integer>>>
    // filterExtractedConnectivitiesByHybridizations(
    // final Map<String, Map<Integer, Map<Integer, Integer>>>
    // extractedConnectivities,
    // final Set<Integer> knownCarbonHybridizations) {
    // // remove hybridization of carbons which we do not expect
    // for (final String atomType : extractedConnectivities.keySet()) {
    // if (atomType.equals("C")) {
    // for (final int hybridization : new
    // HashSet<>(extractedConnectivities.get(atomType)
    // .keySet())) {
    // if (!knownCarbonHybridizations.contains(hybridization)) {
    // extractedConnectivities.get(atomType)
    // .remove(hybridization);
    // }
    // }
    // }
    // }
    //
    // return extractedConnectivities;
    // }
    //
    // @Deprecated
    // public static Map<String, Map<Integer, Set<Integer>>>
    // filterExtractedConnectivitiesByCount(
    // final Map<String, Map<Integer, Map<Integer, Integer>>>
    // extractedConnectivities,
    // final double thresholdElementCount, final boolean onAtomTypeLevel) {
    // final Map<String, Integer> totalCounts =
    // getTotalCounts(extractedConnectivities);
    // final int totalCountsSum = getSum(new HashSet<>(totalCounts.values()));
    // final Map<String, Map<Integer, Set<Integer>>> filteredExtractedConnectivities
    // = new HashMap<>();
    // extractedConnectivities.keySet()
    // .forEach(neighborAtomType -> {
    // int sum = 0;
    // for (final Map.Entry<Integer, Map<Integer, Integer>> entryPerHybridization :
    // extractedConnectivities.get(
    // neighborAtomType)
    // .entrySet()) {
    // for (final Map.Entry<Integer, Integer> entryProtonsCount :
    // extractedConnectivities.get(
    // neighborAtomType)
    // .get(entryPerHybridization.getKey())
    // .entrySet()) {
    // if (onAtomTypeLevel) {
    // sum += entryProtonsCount.getValue();
    // } else if (entryProtonsCount.getValue()
    // / (double) totalCountsSum
    // >= thresholdElementCount) {
    // filteredExtractedConnectivities.putIfAbsent(neighborAtomType,
    // new HashMap<>());
    // filteredExtractedConnectivities.get(neighborAtomType)
    // .putIfAbsent(
    // entryPerHybridization.getKey(),
    // new HashSet<>());
    // filteredExtractedConnectivities.get(neighborAtomType)
    // .get(entryPerHybridization.getKey())
    // .add(entryProtonsCount.getKey());
    // }
    // }
    // }
    // if (onAtomTypeLevel
    // && sum
    // / (double) totalCountsSum
    // >= thresholdElementCount) {
    // filteredExtractedConnectivities.putIfAbsent(neighborAtomType, new
    // HashMap<>());
    // }
    // });
    //
    // return filteredExtractedConnectivities;
    // }
    //
    // @Deprecated
    // private static Map<String, Integer> getTotalCounts(
    // final Map<String, Map<Integer, Map<Integer, Integer>>>
    // extractedConnectivities) {
    // final Map<String, Integer> totalCounts = new HashMap<>();
    // for (final String key1 : extractedConnectivities.keySet()) {
    // totalCounts.putIfAbsent(key1, 0);
    // for (final int key2 : extractedConnectivities.get(key1)
    // .keySet()) {
    // for (final Map.Entry<Integer, Integer> countsEntry :
    // extractedConnectivities.get(key1)
    // .get(key2)
    // .entrySet()) {
    // totalCounts.put(key1, totalCounts.get(key1)
    // + countsEntry.getValue());
    // }
    // }
    // }
    //
    // return totalCounts;
    // }
    //
    // @Deprecated
    // private static int getSum(final Set<Integer> values) {
    // return values.stream()
    // .reduce(0, (total, current) -> total += current);
    // }
}
