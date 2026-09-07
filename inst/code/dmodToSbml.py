#!/usr/bin/env python3
"""Build an SBML Level 3 Version 2 document from a JSON spec produced by
dMod's `exportSbml()`. Mirrors the import helper `sbmlImport.py` but in
the opposite direction. Uses `libsbml.parseL3Formula` to lift the text-form
kinetic laws into MathML ASTs.
"""

import re
import sys
import json
import libsbml


# R writes the natural logarithm as `log`, while an L3 formula reads `log` as
# base 10 and spells the natural one `ln`. Without the rename every `log(x)`
# would come back as `log10(x)`.
_LOG_CALL = re.compile(r"(?<![A-Za-z0-9_.])log\s*\(")


def _l3(formula):
    return _LOG_CALL.sub("ln(", formula)


def _parse(formula, ctx):
    ast = libsbml.parseL3Formula(_l3(formula))
    if ast is None:
        raise ValueError("Could not parse %s: %s -- got %r"
                         % (ctx, libsbml.getLastParseL3Error(), formula))
    return ast


def _check(rc, ctx):
    if rc != libsbml.LIBSBML_OPERATION_SUCCESS:
        raise RuntimeError("libsbml failure at %s: %d" % (ctx, rc))


def build_sbml(spec):
    ns = libsbml.SBMLNamespaces(3, 2)
    document = libsbml.SBMLDocument(ns)
    model = document.createModel()
    _check(model.setId(spec.get("modelId", "dMod_export")), "setId")

    for c in spec.get("compartments", []):
        comp = model.createCompartment()
        _check(comp.setId(c["id"]), "compartment.setId")
        comp.setSize(c.get("size", 1.0))
        comp.setSpatialDimensions(c.get("spatialDimensions", 3))
        comp.setConstant(bool(c.get("constant", True)))
        # A symbolic volume is the compartment's size, not a factor hidden in
        # the kinetic laws: SBML divides by the size, so it has to carry it.
        formula = c.get("sizeAssignment")
        if formula is not None:
            ia = model.createInitialAssignment()
            _check(ia.setSymbol(c["id"]), "initialAssignment.setSymbol")
            ast = libsbml.parseL3Formula(_l3(formula))
            if ast is None:
                raise ValueError(
                    "Could not parse compartment size for %r: %s -- got %r"
                    % (c["id"], libsbml.getLastParseL3Error(), formula))
            ia.setMath(ast)

    for s in spec.get("species", []):
        sp = model.createSpecies()
        _check(sp.setId(s["id"]), "species.setId")
        _check(sp.setCompartment(s["compartment"]), "species.setCompartment")
        # Symbolic initials become InitialAssignments (created below); numeric
        # initials use initialConcentration. SBML lets both coexist -- the
        # InitialAssignment wins at sim time -- but we keep them mutually
        # exclusive for cleanliness on roundtrip.
        amount = bool(s.get("hasOnlySubstanceUnits", False))
        if "initialAssignment" not in s:
            if amount:
                sp.setInitialAmount(s.get("initialAmount", 0.0))
            else:
                sp.setInitialConcentration(s.get("initialConcentration", 0.0))
        sp.setHasOnlySubstanceUnits(amount)
        sp.setBoundaryCondition(False)
        sp.setConstant(False)

    for s in spec.get("species", []):
        formula = s.get("initialAssignment")
        if formula is None:
            continue
        ia = model.createInitialAssignment()
        _check(ia.setSymbol(s["id"]), "initialAssignment.setSymbol")
        ast = libsbml.parseL3Formula(_l3(formula))
        if ast is None:
            raise ValueError(
                "Could not parse initialAssignment for %r: %s -- got %r"
                % (s["id"], libsbml.getLastParseL3Error(), formula))
        ia.setMath(ast)

    for p in spec.get("parameters", []):
        par = model.createParameter()
        _check(par.setId(p["id"]), "parameter.setId")
        par.setValue(p.get("value", 0.0))
        par.setConstant(True)

    for r in spec.get("reactions", []):
        rxn = model.createReaction()
        _check(rxn.setId(r["id"]), "reaction.setId")
        rxn.setReversible(False)
        for reactant in r.get("reactants", []):
            sref = rxn.createReactant()
            sref.setSpecies(reactant["species"])
            sref.setStoichiometry(float(reactant["stoich"]))
            sref.setConstant(True)
        for product in r.get("products", []):
            sref = rxn.createProduct()
            sref.setSpecies(product["species"])
            sref.setStoichiometry(float(product["stoich"]))
            sref.setConstant(True)
        kl = rxn.createKineticLaw()
        ast = libsbml.parseL3Formula(_l3(r["kineticLaw"]))
        if ast is None:
            raise ValueError(
                "Could not parse kinetic law %r: %s"
                % (r["kineticLaw"], libsbml.getLastParseL3Error())
            )
        kl.setMath(ast)

    # --- events ---
    # One <event> per trigger; `time` gives `time >= T`, otherwise the root
    # expression is used as the trigger formula directly.
    for e in spec.get("events", []):
        ev = model.createEvent()
        ev.setUseValuesFromTriggerTime(False)
        trig = ev.createTrigger()
        trig.setInitialValue(False)
        trig.setPersistent(True)
        formula = e.get("trigger")
        ast = libsbml.parseL3Formula(_l3(formula))
        if ast is None:
            raise ValueError(
                "Could not parse event trigger %r: %s"
                % (formula, libsbml.getLastParseL3Error())
            )
        trig.setMath(ast)
        for a in e.get("assignments", []):
            ea = ev.createEventAssignment()
            ea.setVariable(a["variable"])
            m = libsbml.parseL3Formula(_l3(a["formula"]))
            if m is None:
                raise ValueError(
                    "Could not parse event assignment %r: %s"
                    % (a["formula"], libsbml.getLastParseL3Error())
                )
            ea.setMath(m)

    return document


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: %s SPEC.json" % __file__)
        sys.exit(1)

    with open(sys.argv[1]) as fh:
        spec = json.load(fh)

    doc = build_sbml(spec)
    libsbml.writeSBMLToFile(doc, spec["outfile"])
