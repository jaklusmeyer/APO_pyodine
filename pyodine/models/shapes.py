from .base import StaticModel


class LinearStaticModel(StaticModel):
    param_names = ['intercept', 'slope']

    @staticmethod
    def eval(x, params):
        return params['intercept'] + params['slope'] * x

    @staticmethod
    def guess_params(chunk):
        raise NotImplementedError