import pypact as pp

filename = "ve0.out"

with pp.Reader(filename) as output:
#    import pdb; pdb.set_trace()
    for t in output.inventory_data:
        print(t.irradiation_time)
        print(t.currenttime)
        print(t.cooling_time)
