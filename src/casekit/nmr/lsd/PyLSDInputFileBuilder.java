package casekit.nmr.lsd;

import casekit.nmr.lsd.model.ElucidationOptions;
import casekit.nmr.model.nmrdisplayer.Correlation;
import casekit.nmr.model.nmrdisplayer.Data;
import casekit.nmr.model.nmrdisplayer.Link;
import casekit.nmr.utils.Utils;

import java.text.SimpleDateFormat;
import java.util.*;

public class PyLSDInputFileBuilder {

    private static String buildHeader() {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("; PyLSD input file created by webCASE\n");
        final SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd 'at' HH:mm:ss z");
        final Date date = new Date(System.currentTimeMillis());
        stringBuilder.append("; ")
                     .append(formatter.format(date));

        return stringBuilder.toString();
    }

    private static String buildFORM(final String mf, final Map<String, Integer> elementCounts) {
        final StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append("; Molecular Formula: ")
                     .append(mf)
                     .append("\n");
        stringBuilder.append("FORM ");
        elementCounts.forEach((elem, count) -> stringBuilder.append(elem)
                                                            .append(" ")
                                                            .append(count)
                                                            .append(" "));

        return stringBuilder.toString();
    }

    private static String buildPIEC() {
        return "PIEC 1";
    }

    private static String buildELIM(final int elimP1, final int elimP2) {
        return "ELIM "
                + elimP1
                + " "
                + elimP2;
    }

    private static Map<Integer, Object[]> buildIndicesMap(final Data data, final Map<String, Integer> elementCounts) {
        // index in correlation data -> [atom type, indices in PyLSD file...]
        final Map<Integer, Object[]> indicesMap = new HashMap<>();
        // init element indices within correlations with same order as in correlation data input
        final int totalHeavyAtomsCount = elementCounts.entrySet()
                                                      .stream()
                                                      .filter(set -> !set.getKey()
                                                                         .equals("H"))
                                                      .map(Map.Entry::getValue)
                                                      .reduce(0, Integer::sum);
        int heavyAtomIndexInPyLSDFile = 1;
        int protonIndexInPyLSDFile = totalHeavyAtomsCount
                + 1;
        int protonsToInsert;
        Correlation correlation;
        for (int i = 0; i
                < data.getCorrelations()
                      .getValues()
                      .size(); i++) {
            correlation = data.getCorrelations()
                              .getValues()
                              .get(i);
            // set entry for each correlation with consideration of equivalences
            if (correlation.getAtomType()
                           .equals("H")) {
                protonsToInsert = 0;
                for (final Link link : correlation.getLink()) {
                    if (link.getExperimentType()
                            .equals("hsqc")
                            || link.getExperimentType()
                                   .equals("hmqc")) {
                        protonsToInsert += data.getCorrelations()
                                               .getValues()
                                               .get(link.getMatch()
                                                        .get(0))
                                               .getEquivalence();
                    }
                }
                indicesMap.put(i, new Object[1
                        + protonsToInsert]);
                indicesMap.get(i)[0] = correlation.getAtomType();
                for (int j = 0; j
                        < protonsToInsert; j++) {
                    indicesMap.get(i)[1
                            + j] = protonIndexInPyLSDFile;
                    protonIndexInPyLSDFile++;
                }
            } else {
                indicesMap.put(i, new Object[1
                        + correlation.getEquivalence()]);
                indicesMap.get(i)[0] = correlation.getAtomType();
                for (int j = 1; j
                        <= correlation.getEquivalence(); j++) {
                    indicesMap.get(i)[j] = heavyAtomIndexInPyLSDFile;
                    heavyAtomIndexInPyLSDFile++;
                }
            }
        }

        return indicesMap;
    }

    private static String buildMULT(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap,
                                    final Map<Integer, List<Integer>> detectedHybridizations) {
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        List<Integer> hybridizations = new ArrayList<>();
        final StringBuilder hybridizationStringBuilder;
        final StringBuilder attachedProtonsCountStringBuilder;

        if (correlation.getHybridization()
                != null
                && !correlation.getHybridization()
                               .isEmpty()) {
            // if hybridization is already given
            if (correlation.getHybridization()
                           .equals("SP")) {
                hybridizations.add(1);
            } else if (correlation.getHybridization()
                                  .equals("SP2")) {
                hybridizations.add(2);
            } else {
                hybridizations.add(3);
            }
        } else {
            // if hybridization is not given then use the detected ones via MongoDB queries
            if (detectedHybridizations.containsKey(index)) {
                hybridizations = detectedHybridizations.get(index);
            }
        }
        if (hybridizations.isEmpty()) {
            hybridizationStringBuilder = new StringBuilder(
                    Constants.defaultHybridizationMap.get(correlation.getAtomType()));
        } else {
            hybridizationStringBuilder = new StringBuilder();
            if (hybridizations.size()
                    > 1) {
                hybridizationStringBuilder.append("(");
            }
            for (int k = 0; k
                    < hybridizations.size(); k++) {
                hybridizationStringBuilder.append(hybridizations.get(k));
                if (k
                        < hybridizations.size()
                        - 1) {
                    hybridizationStringBuilder.append(" ");
                }
            }
            if (hybridizations.size()
                    > 1) {
                hybridizationStringBuilder.append(")");
            }
        }
        // set attached protons count
        attachedProtonsCountStringBuilder = new StringBuilder();
        // if protons count is given
        if (correlation.getProtonsCount()
                != null
                && !correlation.getProtonsCount()
                               .isEmpty()) {
            if (correlation.getProtonsCount()
                           .size()
                    == 1) {
                attachedProtonsCountStringBuilder.append(correlation.getProtonsCount()
                                                                    .get(0));
            } else {
                attachedProtonsCountStringBuilder.append("(");
                for (final int protonsCount : correlation.getProtonsCount()) {
                    attachedProtonsCountStringBuilder.append(protonsCount)
                                                     .append(" ");
                }
                attachedProtonsCountStringBuilder.deleteCharAt(attachedProtonsCountStringBuilder.length()
                                                                       - 1);
                attachedProtonsCountStringBuilder.append(")");
            }
        } else { // if protons count is not given then set it to default value
            attachedProtonsCountStringBuilder.append(Constants.defaultProtonsCountPerValencyMap.get(
                    Constants.defaultAtomLabelMap.get(correlation.getAtomType())));
        }
        for (int j = 1; j
                < indicesMap.get(index).length; j++) {
            stringBuilder.append("MULT ")
                         .append(indicesMap.get(index)[j])
                         .append(" ")
                         .append(correlation.getAtomType())
                         .append(" ")
                         .append(hybridizationStringBuilder.toString())
                         .append(" ")
                         .append(attachedProtonsCountStringBuilder.toString());
            if (j
                    >= 2) {
                stringBuilder.append("; equivalent to ")
                             .append(indicesMap.get(index)[1]);
            }
            stringBuilder.append("\n");
        }

        return stringBuilder.toString();
    }

    private static String buildHSQC(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hsqc")
                    || link.getExperimentType()
                           .equals("hmqc")) {
                for (final int matchIndex : link.getMatch()) {
                    // for each equivalence of heavy atom and attached equivalent proton
                    for (int k = 1; k
                            < indicesMap.get(index).length; k++) {
                        stringBuilder.append("HSQC ")
                                     .append(indicesMap.get(index)[k])
                                     .append(" ")
                                     .append(indicesMap.get(matchIndex)[k])
                                     .append("\n");
                    }
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildHMBC(final Correlation correlation, final int index, final Data data,
                                    final Map<Integer, Object[]> indicesMap, final int hmbcP3, final int hmbcP4) {
        if (correlation.getAtomType()
                       .equals("H")) {
            return null;
        }
        final String defaultBondDistance = hmbcP3
                + " "
                + hmbcP4;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("hmbc")) {
                for (final int matchIndex : link.getMatch()) {
                    for (int k = 1; k
                            < indicesMap.get(index).length; k++) {
                        for (int l = 1; l
                                < indicesMap.get(matchIndex).length; l++) {
                            // only add an HMBC correlation if there is no direct link via HSQC and the equivalence index is not equal
                            if (!(data.getCorrelations()
                                      .getValues()
                                      .get(matchIndex)
                                      .getAttachment()
                                      .containsKey(correlation.getAtomType())
                                    && data.getCorrelations()
                                           .getValues()
                                           .get(matchIndex)
                                           .getAttachment()
                                           .get(correlation.getAtomType())
                                           .contains(index)
                                    && l
                                    == k)) {
                                uniqueSet.add(indicesMap.get(index)[k]
                                                      + " "
                                                      + indicesMap.get(matchIndex)[l]);
                            }
                        }
                    }
                }
            }
        }

        return uniqueSet.stream()
                        .map(str -> "HMBC "
                                + str
                                + " "
                                + defaultBondDistance
                                + "\n")
                        .reduce("", (strAll, str) -> strAll
                                + str);
    }

    private static String buildCOSY(final Correlation correlation, final int index, final Data data,
                                    final Map<Integer, Object[]> indicesMap, final int cosyP3, final int cosyP4) {
        if (!correlation.getAtomType()
                        .equals("H")) {
            return null;
        }
        final String defaultBondDistance = cosyP3
                + " "
                + cosyP4;
        final Set<String> uniqueSet = new LinkedHashSet<>(); // in case of same content exists multiple times
        for (final Link link : correlation.getLink()) {
            if (link.getExperimentType()
                    .equals("cosy")) {
                for (final int matchIndex : link.getMatch()) {
                    // only add an COSY correlation if the two signals there is not equivalent
                    if (!data.getCorrelations()
                             .getValues()
                             .get(matchIndex)
                             .getId()
                             .equals(correlation.getId())) {
                        for (int k = 1; k
                                < indicesMap.get(index).length; k++) {
                            // only allow COSY values between possible equivalent protons and only one another non-equivalent proton
                            if (indicesMap.get(matchIndex).length
                                    == 2) {
                                uniqueSet.add(indicesMap.get(index)[k]
                                                      + " "
                                                      + indicesMap.get(matchIndex)[1]);
                            }
                        }
                    }
                }
            }
        }

        return uniqueSet.stream()
                        .map(str -> "COSY "
                                + str
                                + " "
                                + defaultBondDistance
                                + "\n")
                        .reduce("", (strAll, str) -> strAll
                                + str);
    }

    private static String buildSHIX(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        if (correlation.getAtomType()
                       .equals("H")
                || correlation.isPseudo()) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        for (int k = 1; k
                < indicesMap.get(index).length; k++) {
            stringBuilder.append("SHIX ")
                         .append(indicesMap.get(index)[k])
                         .append(" ")
                         .append(correlation.getSignal()
                                            .getDelta())
                         .append("\n");
        }

        return stringBuilder.toString();
    }

    private static String buildSHIH(final Correlation correlation, final int index,
                                    final Map<Integer, Object[]> indicesMap) {
        if (!correlation.getAtomType()
                        .equals("H")
                || correlation.isPseudo()) {
            return null;
        }
        final StringBuilder stringBuilder = new StringBuilder();
        // only consider protons which are attached via HSQC/HMQC (pseudo and real links)
        for (final Link link : correlation.getLink()) {
            if ((link.getExperimentType()
                     .equals("hsqc")
                    || link.getExperimentType()
                           .equals("hmqc"))
                    && !link.getMatch()
                            .isEmpty()) { // && !link.isPseudo()
                for (int k = 1; k
                        < indicesMap.get(index).length; k++) {
                    stringBuilder.append("SHIH ")
                                 .append(indicesMap.get(index)[k])
                                 .append(" ")
                                 .append(correlation.getSignal()
                                                    .getDelta())
                                 .append("\n");
                }
            }
        }

        return stringBuilder.toString();
    }

    private static String buildLISTAndPROP(final boolean allowHeteroHeteroBonds) {
        final StringBuilder stringBuilder = new StringBuilder();
        // LIST PROP for hetero hetero bonds allowance
        if (!allowHeteroHeteroBonds) {
            // create hetero atom list automatically
            stringBuilder.append("HETE L1")
                         .append("; list of hetero atoms\n");
            stringBuilder.append("PROP L1 0 L1 -; no hetero-hetero bonds\n");
        }

        return stringBuilder.toString();
    }

    private static String buildFilters(final String[] filterPaths) {
        final StringBuilder stringBuilder = new StringBuilder();
        // DEFF + FEXP -> add filters
        stringBuilder.append("; externally defined filters\n");
        final Map<String, String> filters = new LinkedHashMap<>();
        int counter = 1;
        for (final String filterPath : filterPaths) {
            filters.put("F"
                                + counter, filterPath);
            counter++;
        }

        if (!filters.isEmpty()) {
            filters.forEach((label, filePath) -> stringBuilder.append("DEFF ")
                                                              .append(label)
                                                              .append(" \"")
                                                              .append(filePath)
                                                              .append("\"\n"));
            stringBuilder.append("\n");

            stringBuilder.append("FEXP \"");
            counter = 0;
            for (final String label : filters.keySet()) {
                stringBuilder.append("NOT ")
                             .append(label);
                if (counter
                        < filters.size()
                        - 1) {
                    stringBuilder.append(" and ");
                }
                counter++;
            }
            stringBuilder.append("\"\n");
        }

        return stringBuilder.toString();
    }

    public static String buildPyLSDInputFileContent(final Data data, final String mf,
                                                    final Map<Integer, List<Integer>> detectedHybridizations,
                                                    final ElucidationOptions elucidationOptions) {
        final Map<String, Map<String, Object>> state = data.getCorrelations()
                                                           .getState();
        final boolean hasErrors = state.keySet()
                                       .stream()
                                       .anyMatch(atomType -> state.get(atomType)
                                                                  .containsKey("error")
                                               && !((Map<String, Object>) state.get(atomType)
                                                                               .get("error")).isEmpty());
        if (mf
                != null
                && !hasErrors) {
            final Map<String, Integer> elementCounts = new LinkedHashMap<>(Utils.getMolecularFormulaElementCounts(mf));
            final StringBuilder stringBuilder = new StringBuilder();
            // create header
            stringBuilder.append(buildHeader())
                         .append("\n\n");
            // FORM
            stringBuilder.append(buildFORM(mf, elementCounts))
                         .append("\n\n");
            // PIEC
            stringBuilder.append(buildPIEC())
                         .append("\n\n");
            // ELIM
            if (elucidationOptions.isUseElim()) {
                stringBuilder.append(buildELIM(elucidationOptions.getElimP1(), elucidationOptions.getElimP2()))
                             .append("\n\n");
            }

            final Map<String, List<String>> collection = new LinkedHashMap<>();
            collection.put("MULT", new ArrayList<>());
            collection.put("HSQC", new ArrayList<>());
            collection.put("HMBC", new ArrayList<>());
            collection.put("COSY", new ArrayList<>());
            collection.put("SHIX", new ArrayList<>());
            collection.put("SHIH", new ArrayList<>());
            // index in correlation data -> [atom type, index in PyLSD file]
            final Map<Integer, Object[]> indicesMap = buildIndicesMap(data, elementCounts);

            Correlation correlation;
            for (int i = 0; i
                    < data.getCorrelations()
                          .getValues()
                          .size(); i++) {
                correlation = data.getCorrelations()
                                  .getValues()
                                  .get(i);
                collection.get("MULT")
                          .add(buildMULT(correlation, i, indicesMap, detectedHybridizations));
                collection.get("HSQC")
                          .add(buildHSQC(correlation, i, indicesMap));
                collection.get("HMBC")
                          .add(buildHMBC(correlation, i, data, indicesMap, elucidationOptions.getHmbcP3(),
                                         elucidationOptions.getHmbcP4()));
                collection.get("COSY")
                          .add(buildCOSY(correlation, i, data, indicesMap, elucidationOptions.getCosyP3(),
                                         elucidationOptions.getCosyP4()));
                collection.get("SHIX")
                          .add(buildSHIX(correlation, i, indicesMap));
                collection.get("SHIH")
                          .add(buildSHIH(correlation, i, indicesMap));
            }

            collection.keySet()
                      .forEach(key -> {
                          collection.get(key)
                                    .stream()
                                    .filter(Objects::nonNull)
                                    .forEach(stringBuilder::append);
                          stringBuilder.append("\n");
                      });

            // BOND (interpretation, INADEQUATE, previous assignments) -> input fragments

            // LIST PROP for certain limitations or properties of atoms in lists, e.g. hetero hetero bonds allowance
            stringBuilder.append(buildLISTAndPROP(elucidationOptions.isAllowHeteroHeteroBonds()))
                         .append("\n");
            // DEFF and FEXP as filters (bad lists)
            stringBuilder.append(buildFilters(elucidationOptions.getFilterPaths()))
                         .append("\n");

            return stringBuilder.toString();
        }

        return "";
    }
}
