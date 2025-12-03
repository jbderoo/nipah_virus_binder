python  /path/to/RFDiffusion/run_inference.py \
        inference.output_prefix=outputs \
        inference.input_pdb=cleaned_nipah.pdb \
        'contigmap.contigs=[A1-402/0 0-120/B1-13/0-150/D1-14/5-30/C1-16/0-30]' \
	    'ppi.hotspot_res=[A43,A203]' \
        contigmap.length=150-250 \
	    potentials.guiding_potentials=[\"type:binder_ROG,weight:2,min_dist:15\"] \
	    potentials.guide_scale=2 \
	    potentials.guide_decay='quadratic' \
        inference.num_designs=1 \
        inference.ckpt_override_path=/path/to/models/Complex_base_ckpt.pt

