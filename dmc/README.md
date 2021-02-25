# DMC (Diagram Model Counter)
DMC computes exact literal-weighted model counts of formulas in conjunctive normal form using algebraic decision diagrams
- Developer: Vu Phan

<!-- ####################################################################### -->

## Installation (Linux)

### Prerequisites
- boost 1.66
- cmake 3.11
- g++ 6.4
- gmp 6.2
- make 3.82
- unzip 6.0
- already included in `../common/lib.zip`:
  - [cudd 3.0](https://github.com/ivmai/cudd)
  - [cxxopts 2.2](https://github.com/jarro2783/cxxopts)
  - [sylvan 1.5](https://trolando.github.io/sylvan)

### Command
#### With Singularity 3.5 (slow)
```bash
sudo make
```
#### Without Singularity (quick)
```bash
./SETUP.sh
```

<!-- ####################################################################### -->

## Examples
If you use Singularity, then replace `./dmc` with `singularity run --bind="/:/host" ./dmc.sif` in the following commands.

### Showing command-line options
#### Command
```bash
./dmc
```
#### Output
```
DMC: Diagram Model Counter

Usage:
  dmc [OPTION...]

 Required options:
      --cf arg  cnf file path (input)
      --jf arg  jt file path (to use stdin, type: '--jf -')

 Optional options:
      --wf=arg  weight format in cnf file:
           1    UNWEIGHTED                                        
           2    MINIC2D                                           
           3    CACHET                                            
           4    WCNF                                              
           5    WPCNF                                             Default arg: 5
      --ps=arg  planning strategy:
           f    FIRST_JOIN_TREE                                   Default arg: f
           t    TIMING                                            
      --jw=arg  jt wait duration before jt planner is killed      Default arg: 2.3 (seconds)
      --pf=arg  performance factor (enable with positive float)   Default arg: 0
      --dv=arg  diagram var order heuristic (negate to invert):
           1    APPEARANCE                                        
           2    DECLARATION                                       
           3    RANDOM                                            
           4    MCS                                               Default arg: +4
           5    LEXP                                              
           6    LEXM                                              
           7    MINFILL                                           
      --dp=arg  diagram package:
           1    CUDD                                              Default arg: 1
           2    SYLVAN                                            
      --sw=arg  Sylvan workers (threads)                          Default arg: 0 (auto-detected)
      --jp=arg  join priority:
           a    ARBITRARY                                         
           l    LARGEST_FIRST                                     
           s    SMALLEST_FIRST                                    Default arg: s
      --rs=arg  random seed                                       Default arg: 2020
      --vl=arg  verbosity level:
           0    solution only                                     
           1    parsed info as well                               Default arg: 1
           2    cnf clauses and literal weights as well           
           3    input lines as well                               
```

### Computing model count given cnf file and join tree from jt file
#### Command
```bash
./dmc --cf=../examples/s27_3_2.wpcnf --jf=../examples/s27_3_2.jt.htb --dp=2 --sw=3
```
#### Output
```
c argv: ./dmc --cf=../examples/s27_3_2.wpcnf --jf=../examples/s27_3_2.jt.htb --dp=2 --sw=3

c processing command-line options...
c cnfFilePath                   ../examples/s27_3_2.wpcnf
c jtFilePath                    ../examples/s27_3_2.jt.htb
c weightFormat                  WPCNF
c planningStrategy              FIRST_JOIN_TREE
c jtWaitSeconds                 2.3
c performanceFactor             0
c diagramVarOrder               MCS
c inverseDiagramVarOrder        0
c diagramPackage                SYLVAN
c workerCount                   3
c joinPriority                  SMALLEST_FIRST
c randomSeed                    2020

c processing cnf formula...
c declaredVarCount              20
c apparentVarCount              20
c declaredClauseCount           43
c apparentClauseCount           43

c procressing join tree...
c after 0.001000s, finished processing first join tree (width 10 | ending on/before line 37)
c declaredVarCount              20
c declaredClauseCount           43
c declaredNodeCount             58
c plannerSeconds                0.005

c computing output...
c ------------------------------------------------------------------
s wmc 0.550265
c ------------------------------------------------------------------

c ==================================================================
c seconds                       0.005
c ==================================================================
```

### Computing model count given join tree (from `lg`) and formula (from cnf file)
#### Command
```bash
cnf="../examples/s27_3_2.wpcnf" && ../lg/lg.sif "/solvers/htd-master/bin/htd_main -s 1234567 --opt width --iterations 0 --strategy challenge --print-progress --preprocessing full" < $cnf | ./dmc --cf=$cnf --jf=-
```
#### Output
```
c argv: ./dmc --cf=../examples/s27_3_2.wpcnf --jf=-

c processing command-line options...
c cnfFilePath                   ../examples/s27_3_2.wpcnf
c jtFilePath                    -
c weightFormat                  WPCNF
c planningStrategy              FIRST_JOIN_TREE
c jtWaitSeconds                 2.3
c performanceFactor             0
c diagramVarOrder               MCS
c inverseDiagramVarOrder        0
c diagramPackage                CUDD
c workerCount                   0
c joinPriority                  SMALLEST_FIRST
c randomSeed                    2020

c processing cnf formula...
c declaredVarCount              20
c apparentVarCount              20
c declaredClauseCount           43
c apparentClauseCount           43

c procressing join tree...
c ==================================================================
c getting join tree from stdin... (end input with 'Enter' then 'Ctrl d')
c after 0.094000s, finished processing first join tree (width 20 | ending on/before line 8)
c successfully killed planner process with PID 19525
c ******************************************************************
c MY_WARNING: planner should have been killed; killing it now
c ******************************************************************
c successfully killed planner process with PID 19525
c getting join tree from stdin: done
c ==================================================================
c declaredVarCount              20
c declaredClauseCount           43
c declaredNodeCount             45
c plannerSeconds                0.0202199

c computing output...
c ------------------------------------------------------------------
s wmc 0.550265
c ------------------------------------------------------------------

c ==================================================================
c seconds                       0.095
c ==================================================================
```