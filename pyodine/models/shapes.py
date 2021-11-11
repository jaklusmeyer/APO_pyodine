from .base import StaticModel


class LinearStaticModel(StaticModel):
    param_names = ['intercept', 'slope']

    @staticmethod
    def eval(x, params):
        return params['intercept'] + params['slope'] * x

    @staticmethod
    def guess_params(chunk):
        raise NotImplementedError


class ParabolicStaticModel(StaticModel):
    param_names = ['p0', 'p1', 'p2']
    
    @staticmethod
    def eval(x, params):
        # Is the order of the parameters correct?
        return params['p2'] + params['p1'] * x + params['p0'] * x**2
    
    @staticmethod
    def guess_params(chunk):
        raise NotImplementedError