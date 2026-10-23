from run_cellrank_v1_batch import run_one

failed = [
    ("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Td0_BCG", ['capture9'], 'S1_Td0_BCG'),
    ("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Tm3_BCG", ['capture10', 'capture14'], 'S1_Tm3_BCG'),
    ("SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Td0_BCG", ['capture7'], 'S5_Td0_BCG'),
    ("SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Tm3_BCG", ['capture8'], 'S5_Tm3_BCG'),
]
for uid, caps, label in failed:
    print(f'=== retrying {label} ===')
    run_one(uid, caps, label)
