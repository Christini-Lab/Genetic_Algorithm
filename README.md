# Christini Lab Genetic Algorithm
---

### Introduction
Christini Lab version of the genetic algorithm originally published by
Kumara Sastry at the Illinois Genetic Algorithms Laboratory (IlliGAL).

Reference:
> Single & Multi-Objective Real-Coded Genetic Algorithms Code  
> Author: Kumara Sastry  
> Illinois Genetic Algorithms Laboratory (IlliGAL)  
> Deparment of General Engineering  
> University of Illinois at Urbana-Champaign  
> 104 S. Mathews Ave, Urbana, IL 61801  

This repository is intended to be used as a subtree of a fitting project. Follow
the directions below to add/update the genetic algorithm to your repository.


### Subtree installation:

Enter your local Git project
```sh
cd <PROJECT DIRECTORY>
```

Add remote URL of GA to your local project.
  * **-f** - fetch from remote immediately
  * **GA_remote** - Name of remote, you can change this

```sh
git remote add -f GA_remote https://pbtech-vc.med.cornell.edu/git/christini-lab/Genetic_Algorithm.git
```

Add genetic algorithm as a subtree of the project.
  * **--prefix Genetic_Algorithm/** - prefix denotes the directory you wish to
  put the GA in, you can change this
  * **GA_remote** - Model remote repository set earlier
  * **--squash** - merges all commits into one for cleaner history

```sh
git subtree add --prefix Genetic_Algorithm/ GA_remote --squash
```


### Updating the Genetic_Algorithm subtree:
Fetch any new changes, then pull changes into subtree directory.
  * **--prefix Genetic_Algorithm/** - Specify the model directory
  * **GA_remote** - Model repository you will pull changes from
  * **master** - Pulling in master branch changes, you can specify another
  branch if required

```sh
git subtree pull --prefix Genetic_Algorithm/ GA_remote master --squash
```


### Loading initial populations:
A population file should be formated as follows:
* Number of solutions/individuals to be loaded
* Parameter values (one set of variables for each individual)
* Objective values (for each individual)
* Constraint values (for each individual)
* Penalty values (for each individual)

Example file below contains:
- **2** individuals
- **10** parameter values set to **1.0**
- objective values set to **5**

>2 
>
>1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 5 
>
>1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 5


### Christini Lab modifications:
 * Separated directory structure, which includes examples and documentation
 * Bulk of functions separated out of main function for easier readability and
   development
 * Uses openMP to parallelize evaluations of each individual


### Examples:
 * **sga_nsga** - Original example from IllGAL. Contains 4 different setting files.


### Documentation:
 * **Illigal** - Original report/instructions from IlliGAL.


### Development Notes:

#### Branching:
 * Create new development branch for any unfinished features
 * Features deemed to be necessity will be merged into master branch
 * Optional features will be kept as separate feature branches

#### Code style:
 * Follow 80 column rule
 * No dangling whitespace
 * New code generally follows Google C++ style
