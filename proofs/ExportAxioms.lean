import MathieuVertexOperators
import TadpoleCancellation
open SocrateAI.Moonshine.MathieuVertexOperators
open SocrateAI.StringTheory.TadpoleCancellation

def main : IO Unit := do
  let jsonStr := s!"\{\n  \"h1\": \{\"num\": {h1.num}, \"den\": {h1.den}},\n  \"h2\": \{\"num\": {h2.num}, \"den\": {h2.den}},\n  \"delta12\": \{\"num\": {delta12.num}, \"den\": {delta12.den}},\n  \"delta23\": \{\"num\": {delta23.num}, \"den\": {delta23.den}},\n  \"delta13\": \{\"num\": {delta13.num}, \"den\": {delta13.den}},\n  \"dimA1\": {dimA1},\n  \"dimA2\": {dimA2},\n  \"rNLNum\": {rNLNum},\n  \"rNLDen\": {rNLDen},\n  \"numSupercharges\": {numSupercharges},\n  \"totalD7Charge\": {totalD7Charge},\n  \"totalO7Charge\": {totalO7Charge}\n}"
  
  IO.FS.writeFile "axioms.json" jsonStr
  IO.println "Exported verified axioms to axioms.json"
