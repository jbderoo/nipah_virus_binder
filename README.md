## Introduction

First and foremost, this was a team effort; I would like to not only acknowledge but equally share this effort with both Dr. Christopher D. Snow and (soon to be Dr.) Ashlyn Chen. I had reached out to Chris asking if he had any interest in competing in this competition, and if so would he be willing to donate some computational resources to the endevour. Not only did he agree, he ran several of his own designs that placed in my submitted top 10. A huge thank you to both of them for their time, energy, and scientific insights.

## What is this?
Recently, there was a call for protein binders against the nipah virus G-protein (I will refer to this as just "nipah virus" from here on out) due to its high lethality and recently the WHO listed it as a top-priority for vaccine development. It binds very tightly to the Ephrin B2 and B3 receptors. For more, please see their write up and call here: https://proteinbase.com/competitions/adaptyv-nipah-competition .

This page serves to document my efforts and how I tackled this problem. It also serves as a place to point to when discussing workflows and budding scientists to when talking about how I tackle things. This is an extremely informal white paper on how I approached the problem.

## Structural investigation
It's always worth starting these with both an investigation of the structure, as well as a releveant literature search. The competition gave us a starting structure for us to investigate, 2vsm, a complex of the nipah virus and the human cell surface receptor, Ephrin-B2. ![image of 2vsm](images/2vsm_structure.png)

The binding interaction sits on the pore of the beta (B-) propeller, and reaches up into the pore, creating a strong "double interaction" - a surface to sit on (top half of the B-barrel) and a loop to reach up into the B-propeller and grab onto several charged residues. Looking at the surface representation of the nipah virus lets us easily identify the Ephrin-B2 burying into the pore (__A__), as well as observe the pore itself (__B__). This insight gives us a good idea of the regions to target when we inevitably make our own, competitive binder. ![2vsm pore and interactions](images/2vsm_pore_double_orientation.png)

After a literature search, I identified several structures that might give additional insight into binding of the pore, or adjacent to the pore. Namely, these are structures 7ty0 and 8k0d. 7ty0 has a system of doubled coils (__A__) running along the side of the structure (it also has some sort of dimer binder seated near the pore entrance, but I found nothing of substantial interest in it). These coils are very long, but we can somewhat easily isolate the substructure that is making the most contacts both with the nipah virus and with one another (__B__). ![7ty0](images/double_coils_full.png)

8k0d is a macaca antibody with an extremely long CDR loop that reaches far into the pocket of the nipah virus (__A__). There are a handful of interactions that are very important; here I am highlighting the side chain interactions. We have two tyrosines pinning the bottom of the loop creating strong hydrogen interactions with the nipah virus (__B, C__). While there are other, less strong contacts that occur in the structure, these two coupled with the tyrosine at the very tip of the loop (creating 3 favorable contacts, __D__) make up a majority of the side chain interactions in the loop. There are a handful of other interactions with side chains from the nipah virus interacting with the backbone of the long loop, further reinforcing the interactions made by this loop are favorable enough to splice together into a motif. Finally, a shot of the heavy chain with all 3 previously mentioned residue interactions are shown in (__E__). ![8k0d and subsequent zoom ins](images/macaca_tongue_full_contacts.png)

## Actual computational design
Both BoltzGen (https://github.com/HannesStark/boltzgen) and BindCraft (https://github.com/martinpacesa/BindCraft) have made binder design exceptionally easy. As the new kid on the block, I decided to experiment using BoltzGen the most. I ran 4 different "modes" of BoltzGen: 
  - free design (no motif, no hotspots/binding site)
  - motif grafting (include a motif, no hot spots/binding site)
  - free design with hot spots (no motif, specify residues that are important for binding)
  - motif grafting with hot spots (specify a motif and specific residues that designed binder should interact with).

This was pretty straight forward to do in BoltzGen; for the hot spots I selected a few residues that played a key role in binding to the Ephrin-B2 or the macaca antibody. For example, picking the two arginines that make favorable hyodrgen bonds with the two tyrosines (__B, C__) are great hot spot residues because they are right at the entrance to the pore, and could easily act as sites for the designed binder to interact with. 