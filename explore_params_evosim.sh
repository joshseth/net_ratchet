#!/bin/zsh

# Simulate the evolution of chosen system under a range of parameter schemes. 

# p_mut range:      0.0001 0.001 0.01 0.1
# sigma_mut range:  0.0001 0.001 0.01 0.1
# p_del range:      0.0001 0.001 0.01 0.1
# p_new range:      0.0001 0.001 0.01 0.1

pmut_range=(0.0001 0.001 0.01 0.1)
sigmut_range=(0.0001 0.001 0.01 0.1)
del_range=(0.0001 0.001 0.01 0.1)
add_range=(0.0001 0.001 0.01 0.1)

for pmut in ${pmut_range}
  for sigmut in ${sigmut_range}
    for del in ${del_range}
      if (($del <= $pmut == 1))
        then
          for add in ${add_range}
            if (($add <= pmut == 1)) && (($add <= $del == 1))
            then

              (EVOLDIR=$(./evolve_population.R oscillator 100 10000 $pmut $sigmut $del $add) 
              echo $EVOLDIR

              ./summary_evo.R $EVOLDIR) &
            fi
        fi
