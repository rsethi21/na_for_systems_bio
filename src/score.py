# default approximates mse
def error(predictions, data, substrate_names, obj_calc=None):
    cost = 0
    count = 0
    for substrate_id, substrate_data in data.items():
        for time_point, truth in substrate_data.items():
            prediction = predictions[int(time_point), substrate_names.index(substrate_id)]
            if obj_calc == None:
                cost += (prediction - float(truth))**2
                count += 1
            else:
                cost += obj_calc(prediction, truth)
                count += 1
    return cost/float(count)
