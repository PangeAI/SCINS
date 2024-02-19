
class SCINS:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __repr__(self):
        return self.__dict__.__repr__()

    # Or instead of these different option we can have an argument and a function to return the desired output
    def as_df(self):
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("Pandas is not installed. Please install it to use this function.")
        return pd.DataFrame(self.__dict__, index=[0])

    def as_dict(self):
        return self.__dict__

    def as_list(self):
        return list(self.__dict__.values())

    def as_str(self):
        pass
