from acheron.helpers import model_evaluators

def test_to_float():
    mics = [32, 12, '>32', '>32.0000', '>16']
    floats = [float(i) for i in [32, 12, 32, 32, 16]]

    for i in range(len(mics)):
        assert(model_evaluators.to_float(mics[i]) == floats[i])

def test_is_equiv():
    mics = [32, 16, '>32', '>32.0000', '>16']
    floats = [model_evaluators.to_float(i) for i in [32, 12, 32, 32, 16]]

    for i in range(len(mics)):
        mic_float = model_evaluators.to_float(mics[i])
        if i == 1:
            assert(not model_evaluators.is_equiv(mic_float,floats[i]))
        else:
            assert(model_evaluators.is_equiv(mic_float,floats[i]))

    assert model_evaluators.is_equiv(0.12, 0.125)

def test_to_resistance_type():
    bps = {}
    bps['AMP'] = [8,[16],32]
    bps['AZM'] = [16,[],32]
    bps['CIP'] = [0.06,[0.12,0.25,0.5],1]

    assert model_evaluators.to_resistance_type(1, bps['AMP']) == 'S'
    assert model_evaluators.to_resistance_type(8, bps['AMP']) == 'S'
    assert model_evaluators.to_resistance_type(16, bps['AMP']) == 'I'
    assert model_evaluators.to_resistance_type(32, bps['AMP']) == 'R'
    assert model_evaluators.to_resistance_type(64, bps['AMP']) == 'R'

    assert model_evaluators.to_resistance_type(8, bps['AZM']) == 'S'
    assert model_evaluators.to_resistance_type(16, bps['AZM']) == 'S'
    assert model_evaluators.to_resistance_type(32, bps['AZM']) == 'R'

    assert model_evaluators.to_resistance_type(0.06, bps['CIP']) == 'S'
    assert model_evaluators.to_resistance_type(0.12, bps['CIP']) == 'I'
    assert model_evaluators.to_resistance_type(0.25, bps['CIP']) == 'I'
    assert model_evaluators.to_resistance_type(1, bps['CIP']) == 'R'


def test_find_error_type():
    #find_error_type(predicted, actual, abx)
    assert model_evaluators.find_error_type(1, '32', "AMP") == "Very Major Error"
    assert model_evaluators.find_error_type('1', 32, "AMP") == "Very Major Error"
    assert model_evaluators.find_error_type(16, 32, "AMP") == "Non Error"
    assert model_evaluators.find_error_type(64, 1, "AMP") == "Major Error"

    assert model_evaluators.find_error_type(0.06, 1, "CIP") == "Very Major Error"
    assert model_evaluators.find_error_type(1, 0.06, "CIP") == "Major Error"
    assert model_evaluators.find_error_type(0.12, '0.125', "CIP") =="Non Error"
