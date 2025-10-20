/**
 * ErtlFunctionalGroupsFinder for CDK
 * Copyright (C) 2019 Sebastian Fritsch
 * <p>
 * Source code is available at <https://github.com/zielesny/ErtlFunctionalGroupsFinder>
 * <p>
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * <p>
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package casekit.nmr.fragments.functionalgroup;

import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import org.openscience.cdk.graph.ConnectedComponents;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.GraphUtil.EdgeToBondMap;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;

import java.util.*;

/**
 * Finds and extracts a molecules's functional groups in a purely rule-based
 * manner.
 * <p>
 * This class implements Peter Ertl's algorithm for the automated detection and
 * extraction
 * of functional groups in organic molecules
 * [Ertl P. An algorithm to identify functional groups in organic molecules. J
 * Cheminform. 2017; 9:36.].
 *
 * @author Sebastian Fritsch
 * @version 1.0.0.0
 */
public class ErtlFunctionalGroupsFinder {

    private final static String CARBONYL_C_MARKER = "Carbonyl-C";
    private static final ILoggingTool log = LoggingToolFactory.createLoggingTool(ErtlFunctionalGroupsFinder.class);
    private final Set<Integer> nonmetalAtomicNumbers;
    private final Mode mode;
    private EdgeToBondMap bondMap;
    private int[][] adjList;
    private HashSet<Integer> markedAtoms;
    private HashMap<Integer, Boolean> aromaticHeteroAtoms; // key: atom idx, value: isInGroup
    private Map<IAtom, List<EnvironmentalC>> environmentsMap;

    /**
     * Default constructor for ErtlFunctionalGroupsFinder.
     */
    public ErtlFunctionalGroupsFinder() {
        this(Mode.DEFAULT);
    }

    /**
     * Constructor for ErtlFunctionalGroupsFinder.
     *
     * @param mode working mode (see {@code ErtlFunctionalGroupsFinder.Mode}).
     */
    public ErtlFunctionalGroupsFinder(final Mode mode) {
        this.mode = mode;

        // init non-metal and non-metalloid atom numbers
        this.nonmetalAtomicNumbers = ImmutableSet.of(1, 2, 6, 7, 8, 9, 10, 15, 16, 17, 18, 34, 35, 36, 53, 54, 86);
    }

    private static final boolean isHeteroatom(final IAtom atom) {
        final int atomicNr = atom.getAtomicNumber();
        return atomicNr != 1
                && atomicNr != 6;
    }

    /**
     * Find all functional groups contained in a molecule.
     * <p>
     * NOTE: The input must consist of one connected structure and may not contain
     * charged atoms, metals or metalloids.
     *
     * @param container the molecule which contains the functional groups (may not
     *                  contain charged atoms, metals,
     *                  metalloids or unconnected components!)
     *
     * @return a list with all functional groups found in the molecule.
     */
    public List<IAtomContainer> find(final IAtomContainer container) {
        return this.find(container, true);
    }

    /**
     * Find all functional groups contained in a molecule.
     * <p>
     * NOTE: The input must consist of one connected structure and may not contain
     * charged atoms, metals or metalloids.
     *
     * @param container the molecule which contains the functional groups (may not
     *                  contain charged atoms, metals,
     *                  metalloids or unconnected components!)
     * @param clone     Use 'false' to reuse the input container's bonds and atoms
     *                  in the extraction of the functional
     *                  groups. This may speed up the extraction and lower the
     *                  memory consumption for processing large
     *                  amounts of data but corrupts the original input container.
     *                  Use 'true' to work with a clone and leave the input
     *                  container intact (default).
     *
     * @return a list with all functional groups found in the molecule.
     */
    public List<IAtomContainer> find(final IAtomContainer container, final boolean clone) {
        // work with a clone?
        final IAtomContainer mol;
        if (clone) {
            try {
                mol = container.clone();
            } catch (final CloneNotSupportedException e) {
                throw new IllegalStateException("Atom container could not be cloned");
            }
        } else {
            mol = container;
        }

        // init GraphUtil & EdgeToBondMap
        this.bondMap = EdgeToBondMap.withSpaceFor(mol);
        this.adjList = GraphUtil.toAdjList(mol, this.bondMap);

        this.checkConstraints(mol);

        // atom marking
        this.markAtoms(mol);

        // extract raw groups
        final List<IAtomContainer> groups = this.extractGroups(mol);

        // handle environment
        if (this.mode == Mode.DEFAULT) {
            this.expandGeneralizedEnvironments(groups);
        } else if (this.mode == Mode.NO_GENERALIZATION) {
            this.expandFullEnvironments(groups);
        } else {
            throw new IllegalStateException("Unknown mode.");
        }

        // clear fields
        this.bondMap = null;
        this.adjList = null;
        this.markedAtoms = null;
        this.aromaticHeteroAtoms = null;
        this.environmentsMap = null;

        return groups;
    }

    /**
     * Mark all atoms and store them in a set for further processing.
     *
     * @param molecule Molecule with atoms to mark
     */
    private void markAtoms(final IAtomContainer molecule) {
        if (this.isDbg()) {
            log.debug("########## Starting search for atoms to mark ... ##########");
        }

        // store marked atoms
        this.markedAtoms = Sets.newHashSetWithExpectedSize(molecule.getAtomCount());
        // store aromatic heteroatoms
        this.aromaticHeteroAtoms = new HashMap<>();

        for (int idx = 0; idx < molecule.getAtomCount(); idx++) {
            // skip atoms that already got marked in a previous iteration
            if (this.markedAtoms.contains(idx)) {
                continue;
            }
            final IAtom cAtom = molecule.getAtom(idx);
            // skip aromatic atoms but add them to set
            if (cAtom.isAromatic()) {
                if (isHeteroatom(cAtom)) {
                    this.aromaticHeteroAtoms.put(idx, false);
                }
                continue;
            }

            final int atomicNr = cAtom.getAtomicNumber();

            // if C...
            if (atomicNr == 6) {
                boolean isMarked = false; // to detect if foor loop ran with or without marking the C atom
                int oNSCounter = 0; // count for the number of connected O, N & S atoms
                for (final int connectedIdx : this.adjList[idx]) {
                    final IAtom connectedAtom = molecule.getAtom(connectedIdx);
                    final IBond connectedBond = this.bondMap.get(idx, connectedIdx);

                    // if connected to Heteroatom or C in aliphatic double or triple bond...
                    // [CONDITIONS 2.1 & 2.2]
                    if (connectedAtom.getAtomicNumber() != 1
                            && ((connectedBond.getOrder() == Order.DOUBLE
                                    || connectedBond.getOrder() == Order.TRIPLE)
                                    && !connectedBond.isAromatic())) {

                        // set the connected atom as marked
                        if (this.markedAtoms.add(connectedIdx)) {
                            final String connectedAtomCondition = connectedAtom.getAtomicNumber() == 6
                                    ? "2.1/2.2"
                                    : "1";
                            if (this.isDbg()) {
                                log.debug(String.format("Marking Atom #%d (%s) - Met condition %s", connectedIdx,
                                        connectedAtom.getSymbol(), connectedAtomCondition));
                            }
                        }

                        // set the current atom as marked and break out of connected atoms
                        if (this.isDbg()) {
                            log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.1/2.2", idx,
                                    cAtom.getSymbol()));
                        }
                        isMarked = true;

                        // but check for carbonyl-C before break
                        if (connectedAtom.getAtomicNumber() == 8
                                && connectedBond.getOrder() == Order.DOUBLE
                                && this.adjList[idx].length == 3) {
                            if (this.isDbg()) {
                                log.debug("                     - was flagged as Carbonly-C");
                            }
                            cAtom.setProperty(CARBONYL_C_MARKER, true);
                        }

                        break;
                    }
                    // if connected to O/N/S in single bond...
                    else if ((connectedAtom.getAtomicNumber() == 7
                            || connectedAtom.getAtomicNumber() == 8
                            || connectedAtom.getAtomicNumber() == 16)
                            && connectedBond.getOrder() == Order.SINGLE) {
                        // if connected O/N/S is not aromatic...
                        if (!connectedAtom.isAromatic()) {
                            // set the connected O/N/S atom as marked
                            if (this.isDbg()) {
                                log.debug(String.format("Marking Atom #%d (%s) - Met condition 1", connectedIdx,
                                        connectedAtom.getSymbol()));
                            }
                            this.markedAtoms.add(connectedIdx);

                            // if "acetal C" (2+ O/N/S in single bonds connected to sp3-C)... [CONDITION
                            // 2.3]
                            boolean isAllSingleBonds = true;
                            for (final int connectedInSphere2Idx : this.adjList[connectedIdx]) {
                                final IBond sphere2Bond = this.bondMap.get(connectedIdx, connectedInSphere2Idx);
                                if (sphere2Bond.getOrder() != Order.SINGLE) {
                                    isAllSingleBonds = false;
                                    break;
                                }
                            }
                            if (isAllSingleBonds) {
                                oNSCounter++;
                                if (oNSCounter > 1
                                        && this.adjList[idx].length
                                                + cAtom.getImplicitHydrogenCount() == 4) {
                                    // set as marked and break out of connected atoms
                                    if (this.isDbg()) {
                                        log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.3", idx,
                                                cAtom.getSymbol()));
                                    }
                                    isMarked = true;
                                    break;
                                }
                            }
                        }
                        // if part of oxirane, aziridine and thiirane ring... [CONDITION 2.4]
                        for (final int connectedInSphere2Idx : this.adjList[connectedIdx]) {
                            final IAtom connectedInSphere2Atom = molecule.getAtom(connectedInSphere2Idx);
                            if (connectedInSphere2Atom.getAtomicNumber() == 6) {
                                for (final int connectedInSphere3Idx : this.adjList[connectedInSphere2Idx]) {
                                    final IAtom connectedInSphere3Atom = molecule.getAtom(connectedInSphere3Idx);
                                    if (connectedInSphere3Atom.equals(cAtom)) {
                                        // set connected atoms as marked
                                        if (this.isDbg()) {
                                            log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4",
                                                    connectedInSphere2Idx,
                                                    connectedInSphere2Atom.getSymbol()));
                                        }
                                        if (this.isDbg()) {
                                            log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4",
                                                    connectedInSphere3Idx,
                                                    connectedInSphere3Atom.getSymbol()));
                                        }
                                        this.markedAtoms.add(connectedInSphere2Idx);
                                        this.markedAtoms.add(connectedInSphere3Idx);
                                        // set current atom as marked and break out of connected atoms
                                        if (this.isDbg()) {
                                            log.debug(String.format("Marking Atom #%d (%s) - Met condition 2.4", idx,
                                                    cAtom.getSymbol()));
                                        }
                                        isMarked = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                if (isMarked) {
                    this.markedAtoms.add(idx);
                    continue;
                }
                // if none of the conditions 2.X apply, we have an unmarked C (not relevant
                // here)
            }
            // if H...
            else if (atomicNr == 1) {
                // convert to implicit H
                final IAtom connectedAtom;
                try {
                    connectedAtom = molecule.getAtom(this.adjList[idx][0]);
                } catch (final ArrayIndexOutOfBoundsException e) {
                    break;
                }

                if (connectedAtom.getImplicitHydrogenCount() == null) {
                    connectedAtom.setImplicitHydrogenCount(1);
                } else {
                    connectedAtom.setImplicitHydrogenCount(connectedAtom.getImplicitHydrogenCount()
                            + 1);
                }
                continue;
            }
            // if heteroatom... (CONDITION 1)
            else {
                if (this.isDbg()) {
                    log.debug(String.format("Marking Atom #%d (%s) - Met condition 1", idx, cAtom.getSymbol()));
                }
                this.markedAtoms.add(idx);
                continue;
            }
        }
        if (this.isDbg()) {
            log.debug(String.format("########## End of search. Marked %d/%d atoms. ##########", this.markedAtoms.size(),
                    molecule.getAtomCount()));
        }
    }

    /**
     * Searches the molecule for groups of connected marked atoms and extracts each
     * as a new functional group.
     * The extraction process includes marked atom's "environments". Connected H's
     * are captured implicitly.
     *
     * @param molecule the molecule which contains the functional groups
     *
     * @return a list of all functional groups (including "environments") extracted
     *         from the molecule
     */
    private List<IAtomContainer> extractGroups(final IAtomContainer molecule) {
        if (this.isDbg()) {
            log.debug("########## Starting identification & extraction of functional groups... ##########");
        }

        this.environmentsMap = Maps.newHashMapWithExpectedSize(molecule.getAtomCount());
        final int[] atomIdxToFGMap = new int[molecule.getAtomCount()];
        Arrays.fill(atomIdxToFGMap, -1);
        int fGroupIdx = -1;

        while (!this.markedAtoms.isEmpty()) {
            // search for another functional group
            fGroupIdx++;

            // get next markedAtom as the starting node for the search
            final int beginIdx = this.markedAtoms.iterator()
                    .next();
            if (this.isDbg()) {
                log.debug(String.format("Searching new functional group from atom #%d (%s)...", beginIdx,
                        molecule.getAtom(beginIdx)
                                .getSymbol()));
            }

            // do a BFS from there
            final Queue<Integer> queue = new ArrayDeque<>();
            queue.add(beginIdx);

            while (!queue.isEmpty()) {
                final int currentIdx = queue.poll();

                // we are only interested in marked atoms that are not yet included in a group
                if (!this.markedAtoms.contains(currentIdx)) {
                    continue;
                }

                // if it isn't...
                final IAtom currentAtom = molecule.getAtom(currentIdx);
                if (this.isDbg()) {
                    log.debug(String.format("	visiting marked atom: #%d (%s)", currentIdx, currentAtom.getSymbol()));
                }

                // add its index to the functional group
                atomIdxToFGMap[currentIdx] = fGroupIdx;
                // also scratch the index from markedAtoms
                this.markedAtoms.remove(currentIdx);

                // and take look at the connected atoms
                final List<EnvironmentalC> currentEnvironment = new ArrayList<>();
                for (final int connectedIdx : this.adjList[currentIdx]) {
                    // add connected marked atoms to queue
                    if (this.markedAtoms.contains(connectedIdx)) {
                        queue.add(connectedIdx);
                        continue;
                    }

                    // ignore already handled connected atoms
                    if (atomIdxToFGMap[connectedIdx] >= 0) {
                        continue;
                    }

                    // add unmarked connected aromatic heteroatoms
                    final IAtom connectedAtom = molecule.getAtom(connectedIdx);
                    if (isHeteroatom(connectedAtom)
                            && connectedAtom.isAromatic()) {
                        if (this.isDbg()) {
                            log.debug("	   added connected aromatic heteroatom "
                                    + connectedAtom.getSymbol());
                        }
                        atomIdxToFGMap[connectedIdx] = fGroupIdx;
                        // note that this aromatic heteroatom has been added to a group
                        this.aromaticHeteroAtoms.put(connectedIdx, true);
                    }

                    // add unmarked connected atoms to current marked atom's environment
                    final IBond connectedBond = this.bondMap.get(currentIdx, connectedIdx);

                    final EnvironmentCalCType type;
                    if (connectedAtom.getAtomicNumber() == 6) {
                        if (connectedAtom.isAromatic()) {
                            type = EnvironmentCalCType.C_AROMATIC;
                        } else {
                            type = EnvironmentCalCType.C_ALIPHATIC;
                        }
                    } else {
                        // aromatic heteroatom, so just ignore
                        continue;
                    }
                    currentEnvironment.add(new EnvironmentalC(type, connectedBond,
                            connectedBond.getBegin() == connectedAtom
                                    ? 0
                                    : 1));
                }
                this.environmentsMap.put(currentAtom, currentEnvironment);

                // debug logging
                if (this.isDbg()) {
                    int cAromCount = 0, cAliphCount = 0;
                    for (final EnvironmentalC comp : currentEnvironment) {
                        if (comp.getType() == EnvironmentCalCType.C_AROMATIC) {
                            cAromCount++;
                        } else if (comp.getType() == EnvironmentCalCType.C_ALIPHATIC) {
                            cAliphCount++;
                        }
                    }
                    log.debug(String.format(
                            "	   logged marked atom's environment: C_ar:%d, C_al:%d (and %d implicit hydrogens)",
                            cAromCount, cAliphCount, currentAtom.getImplicitHydrogenCount()));
                }
            }

            if (this.isDbg()) {
                log.debug("	search completed.");
            }
        }

        // also create FG for lone aromatic heteroatoms, not connected to a FG yet.
        for (final int atomIdx : this.aromaticHeteroAtoms.keySet()) {
            if (!this.aromaticHeteroAtoms.get(atomIdx)) {
                fGroupIdx++;
                atomIdxToFGMap[atomIdx] = fGroupIdx;
                if (this.isDbg()) {
                    log.debug("Created FG for lone aromatic heteroatom: "
                            + molecule.getAtom(atomIdx)
                                    .getSymbol());
                }
            }
        }

        final List<IAtomContainer> fGs = this.partitionIntoGroups(molecule, atomIdxToFGMap, fGroupIdx
                + 1);

        if (this.isDbg()) {
            log.debug(String.format("########## Found & extracted %d functional groups. ##########", fGroupIdx
                    + 1));
        }
        return fGs;
    }

    /**
     * Generalizes the full environments of functional groups, providing a good
     * balance between preserving
     * meaningful detail and generalization.
     *
     * @param fGroups the list of functional groups including "environments"
     */
    private void expandGeneralizedEnvironments(final List<IAtomContainer> fGroups) {
        if (this.isDbg()) {
            log.debug("########## Starting generalization of functional groups... ##########");
        }

        for (final IAtomContainer fGroup : fGroups) {
            final int atomCount = fGroup.getAtomCount();

            if (this.isDbg()) {
                log.debug(String.format("Generalizing functional group (%d atoms)...", atomCount));
            }

            // prechecking for special cases...
            if (fGroup.getAtomCount() == 1) {
                final IAtom atom = fGroup.getAtom(0);
                final List<EnvironmentalC> environment = this.environmentsMap.get(atom);

                if (environment != null) {
                    final int envCCount = environment.size();

                    // for H2N-C_env & HO-C_env -> do not replace H & C_env by R!
                    if ((atom.getAtomicNumber() == 8
                            && envCCount == 1)
                            || (atom.getAtomicNumber() == 7
                                    && envCCount == 1)) {
                        if (this.isDbg()) {
                            log.debug(String.format(
                                    "   - found single atomic N or O FG with one env. C. Expanding environment...",
                                    atom.getSymbol()));
                        }
                        this.expandEnvironment(atom, fGroup);

                        final int hCount = atom.getImplicitHydrogenCount();
                        if (hCount != 0) {
                            if (this.isDbg()) {
                                log.debug(String.format("   - adding %d hydrogens...", hCount));
                            }
                            this.addHydrogens(atom, hCount, fGroup);
                            atom.setImplicitHydrogenCount(0);
                        }
                        continue;
                    }
                    // for HN-(C_env)-C_env & HS-C_env -> do not replace H by R! (only C_env!)
                    if ((atom.getAtomicNumber() == 7
                            && envCCount == 2)
                            || (atom.getAtomicNumber() == 16
                                    && envCCount == 1)) {
                        if (this.isDbg()) {
                            log.debug("   - found sec. amine or simple thiol");
                        }
                        final int hCount = atom.getImplicitHydrogenCount();
                        if (hCount != 0) {
                            if (this.isDbg()) {
                                log.debug(String.format("   - adding %d hydrogens...", hCount));
                            }
                            this.addHydrogens(atom, hCount, fGroup);
                            atom.setImplicitHydrogenCount(0);
                        }
                        if (this.isDbg()) {
                            log.debug("   - expanding environment...");
                        }
                        this.expandEnvironmentGeneralized(atom, fGroup);
                        continue;
                    }
                } else if (isHeteroatom(atom)) {
                    final int rAtomCount = atom.getValency();
                    final Integer hCount = atom.getImplicitHydrogenCount();
                    if (hCount != null
                            && hCount != 0) {
                        atom.setImplicitHydrogenCount(0);
                    }
                    final String atomTypeName = atom.getAtomTypeName();
                    if (this.isDbg()) {
                        log.debug(String.format(
                                "   - found single aromatic heteroatom (%s, Atomtype %s). Adding %d R-Atoms...",
                                atom.getSymbol(), atomTypeName, rAtomCount));
                    }
                    this.addRAtoms(atom, rAtomCount, fGroup);
                    continue;
                }
            }

            // get atoms to process
            final List<IAtom> fGroupAtoms = Lists.newArrayList(fGroup.atoms());

            // process atoms...
            for (final IAtom atom : fGroupAtoms) {
                final List<EnvironmentalC> environment = this.environmentsMap.get(atom);

                if (environment == null) {
                    if (atom.getImplicitHydrogenCount() != 0) {
                        atom.setImplicitHydrogenCount(0);
                    }
                    final int rAtomCount = atom.getValency()
                            - 1;
                    if (this.isDbg()) {
                        log.debug(String.format("   - found connected aromatic heteroatom (%s). Adding %d R-Atoms...",
                                atom.getSymbol(), rAtomCount));
                    }
                    this.addRAtoms(atom, rAtomCount, fGroup);
                }

                // processing carbons...
                if (atom.getAtomicNumber() == 6) {
                    if (atom.getProperty(CARBONYL_C_MARKER) == null) {
                        if (atom.getImplicitHydrogenCount() != 0) {
                            atom.setImplicitHydrogenCount(0);
                        }
                        if (this.isDbg()) {
                            log.debug("   - ignoring environment for marked carbon atom");
                        }
                        continue;
                    } else {
                        if (this.isDbg()) {
                            log.debug("   - found carbonyl-carbon. Expanding environment...");
                        }
                        this.expandEnvironmentGeneralized(atom, fGroup);
                        continue;
                    }
                }
                // processing heteroatoms...
                else {
                    if (this.isDbg()) {
                        log.debug(String.format("   - found heteroatom (%s). Expanding environment...",
                                atom.getSymbol()));
                    }
                    this.expandEnvironmentGeneralized(atom, fGroup);
                    continue;
                }
            }
        }

        if (this.isDbg()) {
            log.debug("########## Generalization of functional groups completed. ##########");
        }
    }

    /**
     * Expands the full environments of functional groups, converted into atoms and
     * bonds.
     *
     * @param fGroups the list of functional groups including "environments"
     */
    private void expandFullEnvironments(final List<IAtomContainer> fGroups) {
        if (this.isDbg()) {
            log.debug("########## Starting expansion of full environments for functional groups... ##########");
        }

        for (final IAtomContainer fGroup : fGroups) {
            final int atomCount = fGroup.getAtomCount();
            if (this.isDbg()) {
                log.debug(String.format("Expanding environment on functional group (%d atoms)...", atomCount));
            }

            for (int i = 0; i < atomCount; i++) {
                final IAtom atom = fGroup.getAtom(i);

                if (this.isDbg()) {
                    log.debug(String.format(" - Atom #%d:%   - Expanding environment...", i));
                }
                this.expandEnvironment(atom, fGroup);

                final int hCount = atom.getImplicitHydrogenCount();
                if (hCount != 0) {
                    if (this.isDbg()) {
                        log.debug(String.format("   - adding %d hydrogens...", hCount));
                    }
                    this.addHydrogens(atom, hCount, fGroup);
                    atom.setImplicitHydrogenCount(0);
                }
            }
        }

        if (this.isDbg()) {
            log.debug("########## Expansion of full environments for functional groups completed. ##########");
        }
    }

    private void expandEnvironment(final IAtom atom, final IAtomContainer container) {
        final List<EnvironmentalC> environment = this.environmentsMap.get(atom);

        if (environment == null
                || environment.isEmpty()) {
            if (this.isDbg()) {
                log.debug("		found no environment to expand.");
            }
            return;
        }

        int cAromCount = 0, cAliphCount = 0;
        for (final EnvironmentalC envC : environment) {
            final IAtom cAtom = atom.getBuilder()
                    .newInstance(IAtom.class, "C");
            cAtom.setAtomTypeName("C");
            cAtom.setImplicitHydrogenCount(0);
            if (envC.getType() == EnvironmentCalCType.C_AROMATIC) {
                cAtom.setIsAromatic(true);
                cAromCount++;
            } else {
                cAliphCount++;
            }

            final IBond bond = envC.createBond(atom, cAtom);

            container.addAtom(cAtom);
            container.addBond(bond);
        }

        if (this.isDbg()) {
            log.debug(String.format("		expanded environment: %dx C_ar and %dx C_al", cAromCount, cAliphCount));
        }
    }

    // only call this on marked heteroatoms / carbonyl-C's!
    private void expandEnvironmentGeneralized(final IAtom atom, final IAtomContainer container) {

        final List<EnvironmentalC> environment = this.environmentsMap.get(atom);

        if (environment == null) {
            if (this.isDbg()) {
                log.debug("		found no environment to expand.");
            }
            return;
        }

        int rAtomCount = environment.size();
        final int rAtomsForCCount = rAtomCount;
        if (atom.getAtomicNumber() == 8
                && atom.getImplicitHydrogenCount() == 1) {
            this.addHydrogens(atom, 1, container);
            atom.setImplicitHydrogenCount(0);
            if (this.isDbg()) {
                log.debug("		expanded hydrogen on connected OH-Group");
            }
        } else if (isHeteroatom(atom)) {
            rAtomCount += atom.getImplicitHydrogenCount();
        }
        this.addRAtoms(atom, rAtomCount, container);

        if (atom.getImplicitHydrogenCount() != 0) {
            atom.setImplicitHydrogenCount(0);
        }

        if (this.isDbg()) {
            log.debug(String.format("		expanded environment: %dx R-atom (incl. %d for H replacement)", rAtomCount,
                    rAtomCount
                            - rAtomsForCCount));
        }
    }

    private final boolean isNonmetal(final IAtom atom) {
        return this.nonmetalAtomicNumbers.contains(atom.getAtomicNumber());
    }

    private void addHydrogens(final IAtom atom, final int number, final IAtomContainer container) {
        for (int i = 0; i < number; i++) {
            final IAtom hydrogen = atom.getBuilder()
                    .newInstance(IAtom.class, "H");
            hydrogen.setAtomTypeName("H");
            hydrogen.setImplicitHydrogenCount(0);

            container.addAtom(hydrogen);
            container.addBond(atom.getBuilder()
                    .newInstance(IBond.class, atom, hydrogen, Order.SINGLE));
        }
    }

    private void addRAtoms(final IAtom atom, final int number, final IAtomContainer container) {
        for (int i = 0; i < number; i++) {
            final IPseudoAtom rAtom = atom.getBuilder()
                    .newInstance(IPseudoAtom.class, "R");
            rAtom.setAttachPointNum(1);
            rAtom.setImplicitHydrogenCount(0);

            container.addAtom(rAtom);
            container.addBond(atom.getBuilder()
                    .newInstance(IBond.class, atom, rAtom, Order.SINGLE));
        }
    }

    private List<IAtomContainer> partitionIntoGroups(final IAtomContainer sourceContainer, final int[] atomIdxToFGMap,
            final int fGroupCount) {
        final List<IAtomContainer> groups = new ArrayList<>(fGroupCount);
        for (int i = 0; i < fGroupCount; i++) {
            groups.add(sourceContainer.getBuilder()
                    .newInstance(IAtomContainer.class));
        }

        final Map<IAtom, IAtomContainer> atomtoFGMap = Maps.newHashMapWithExpectedSize(sourceContainer.getAtomCount());

        // atoms
        for (int atomIdx = 0; atomIdx < sourceContainer.getAtomCount(); atomIdx++) {
            final int fGroupId = atomIdxToFGMap[atomIdx];

            if (fGroupId == -1) {
                continue;
            }

            final IAtom atom = sourceContainer.getAtom(atomIdx);
            final IAtomContainer myGroup = groups.get(fGroupId);
            myGroup.addAtom(atom);
            atomtoFGMap.put(atom, myGroup);
        }

        // bonds
        for (final IBond bond : sourceContainer.bonds()) {
            final IAtomContainer beginGroup = atomtoFGMap.get(bond.getBegin());
            final IAtomContainer endGroup = atomtoFGMap.get(bond.getEnd());

            if (beginGroup == null
                    || endGroup == null
                    || beginGroup != endGroup) {
                continue;
            }

            beginGroup.addBond(bond);
        }

        // single electrons
        for (final ISingleElectron electron : sourceContainer.singleElectrons()) {
            final IAtomContainer group = atomtoFGMap.get(electron.getAtom());
            if (group != null) {
                group.addSingleElectron(electron);
            }
        }

        // lone pairs
        for (final ILonePair lonePair : sourceContainer.lonePairs()) {
            final IAtomContainer group = atomtoFGMap.get(lonePair.getAtom());
            if (group != null) {
                group.addLonePair(lonePair);
            }
        }

        return groups;
    }

    private boolean isDbg() {
        return log.isDebugEnabled();
    }

    private boolean checkConstraints(final IAtomContainer molecule) {
        for (final IAtom atom : molecule.atoms()) {
            if (atom.getFormalCharge() != null
                    && atom.getFormalCharge() != 0) {
                throw new IllegalArgumentException("Input molecule must not contain any charges.");
            }
            if (!this.isNonmetal(atom)) {
                throw new IllegalArgumentException("Input molecule must not contain metals or metalloids.");
            }
            if (atom.getImplicitHydrogenCount() == null) {
                atom.setImplicitHydrogenCount(0);
            }
        }

        final ConnectedComponents cc = new ConnectedComponents(this.adjList);
        if (cc.nComponents() != 1) {
            throw new IllegalArgumentException("Input molecule must consist of only a single connected stucture.");
        }

        return true;
    }

    /**
     * Defines the working mode.
     */
    public enum Mode {
        /**
         * Default mode including the generalization step.
         */
        DEFAULT,
        /**
         * Skips the generalization step. Functional groups will keep their full
         * "environment".
         */
        NO_GENERALIZATION
    }

    private enum EnvironmentCalCType {
        C_AROMATIC, C_ALIPHATIC
    }

    /**
     * Describes one carbon atom in the environment of a marked atom. It can either
     * be aromatic
     * or aliphatic and also contains a clone of its connecting bond.
     */
    private class EnvironmentalC {
        private final EnvironmentCalCType type;
        private final int bondIndex;
        private final Order bondOrder;
        private final IBond.Stereo bondStereo;
        private final boolean[] bondFlags;

        public EnvironmentalC(final EnvironmentCalCType type, final IBond bond, final int indexInBond) {
            this.type = type;

            this.bondIndex = indexInBond;
            this.bondOrder = bond.getOrder();
            this.bondStereo = bond.getStereo();
            this.bondFlags = bond.getFlags();
        }

        public EnvironmentCalCType getType() {
            return this.type;
        }

        public IBond createBond(final IAtom targetAtom, final IAtom cAtom) {
            final IBond bond = targetAtom.getBuilder()
                    .newInstance(IBond.class);
            if (this.bondIndex == 0) {
                bond.setAtoms(new IAtom[] { cAtom, targetAtom });
            } else {
                bond.setAtoms(new IAtom[] { targetAtom, cAtom });
            }
            bond.setOrder(this.bondOrder);
            bond.setStereo(this.bondStereo);
            bond.setFlags(this.bondFlags);

            return bond;
        }
    }
}