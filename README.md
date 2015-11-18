Christini Lab version of the genetic algorithm originally published by
Kumara Sastry at the Illinois Genetic Algorithms Laboratory (IlliGAL).

Reference:  
  Single & Multi-Objective Real-Coded Genetic Algorithms Code  
  Author: Kumara Sastry  
  Illinois Genetic Algorithms Laboratory (IlliGAL)  
  Deparment of General Engineering  
  University of Illinois at Urbana-Champaign  
  104 S. Mathews Ave, Urbana, IL 61801  

This repository is intended to be used as a subtree of a fitting project. Follow
the directions below to add/update the genetic algorithm to your repository.
Subtree installation:  
  # Enter your local Git project
```
cd <PROJECT DIRECTORY>
```

  # Add remote URL of GA to your local project
```
git remote add -f Genetic_Algorithm https://pbtech-vc.med.cornell.edu/git/christini-lab/Genetic_Algorithm.git
```
  # Merge master branch of Genetic Algorithm into local project. This will not
    change any of your files locally. You can replace master with any other
    branch you wish to use.
```
git merge -s ours --no-commit Genetic_Algorithm/master
```

  # Create a new directory called "Genetic_Algorithm" (You can change this to
    whatever you wish after --prefix=) and copy Git history of master branch of
    the Genetic_Algorithm repository.
```
  git read-tree --prefix=Genetic_Algorithm/ -u Genetic_Algorithm/master
```

  # Commit the changes
```
  git commit -m "Genetic_Algorithm merged as subtree"
```

Updating the Genetic_Algorithm subtree:
  # Pull changes from Genetic_Algorithm master branch. If using another branch,
    change master accordingly.
```
  git pull -s subtree Genetic_Algorithm master
```

Christini Lab modifications:  
  * Separated directory structure, which includes examples and documentation
  * Bulk of functions separated out of main function for easier readability and
    development
  * Uses openMP to parallelize evaluations of each individual

Examples:  
  sga_nsga: Original example from IllGAL. Contains 4 different setting files.

Documentation:  
  Illigal: Original report/instructions from IlliGAL.

Development Notes:  
  Branching:  
    * Create new development branch for any unfinished features
    * Features deemed to be necessity will be merged into master branch
    * Optional features will be kept as separate feature branches
  Code style:  
    * Follow 80 column rule
    * No dangling whitespace
    * New code generally follows Google C++ style
