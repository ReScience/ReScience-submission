import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Flatten, Conv1D, MaxPooling1D, Reshape, UpSampling1D
from tensorflow.keras.models import Model
from tensorflow.keras import losses
from sklearn.base import BaseEstimator


class autoenconder(BaseEstimator):
    """ Template
    TO-DO: detailed explanation

    Parameters
    ----------




    """ 
    def __init__(self, epochs=10, batch_size=32,
                 name_dataset=None, value_encoding_dim = 2, 
                 type_loss = 'l1'):
        """ Template
        TO-DO: detailed explanation
        
        Parameters
        ----------


        

        """        
        # auto-enconder parameters
        self.value_encoding_dim = value_encoding_dim
        self.batch_size = batch_size
        self.epochs = epochs
        self.type_loss = type_loss
        
        

        # information about when is use
        self.name_dataset = name_dataset

        # information about history train
        self.history = []

        # saving method
        self.method_autoenconder = []
        self.method_enconder = []
        
        self.build_auto_enconder()
        
    def build_auto_enconder(self):
        """ Template
        TO-DO: detailed explanation
        
        Parameters
        ----------

        value_encoding_dim : TO-DO
        type_loss: TO-DO
        

        """        
        if(self.type_loss == 'l1'):
            fun_loss = losses.mean_absolute_error
        else:
            raise ValueError("Loss function not yet implemented.")
        
        
        original_signal = Input(shape=(4096, 1))
        
        enconded = Conv1D(kernel_size=3, filters=16,
                          padding='same', activation='relu')(original_signal)
        
        enconded = MaxPooling1D(pool_size=2)(enconded)
        
        enconded = Conv1D(kernel_size=3, filters=32,
                          padding='same', activation='relu')(enconded)
        
        enconded = MaxPooling1D(pool_size=2)(enconded)
        
        enconded = Conv1D(kernel_size=3, filters=64,
                          padding='same', activation='relu')(enconded)
        
        enconded = MaxPooling1D(pool_size=2)(enconded)
        
        enconded = Flatten()(enconded)

        enconded = Dense(self.value_encoding_dim, activation='relu')(enconded)

        decoded = Dense(512*64, activation='relu', use_bias=False)(enconded)

        decoded = Reshape((512, 64))(decoded)

        decoded = Conv1D(kernel_size=3, filters=64,
                         padding='same', activation='relu')(decoded)
        decoded = UpSampling1D()(decoded)
        
        decoded = Conv1D(kernel_size=3, filters=32,
                         padding='same', activation='relu')(decoded)
        
        decoded = UpSampling1D()(decoded)
        
        decoded = Conv1D(kernel_size=3, filters=16,
                         padding='same', activation='relu')(decoded)
        decoded = UpSampling1D()(decoded)
        
        decoded = Conv1D(kernel_size=3, filters=1,
                         padding='same', activation='sigmoid')(decoded)
        
        encoder = Model(original_signal, enconded, name='encoder')

        autoencoder = Model(original_signal, decoded, 
                            name='autoenconder_'+str(self.value_encoding_dim))
        
        autoencoder.compile(optimizer='adam', loss=fun_loss, metrics=['accuracy'])
        
        self.method_autoenconder = autoencoder
        self.method_enconder = encoder

    def fit(self, X_train, X_validation):
        """ Template
        TO-DO: detailed explanation
        
        Parameters
        ----------


        

        """ 

        # Training auto-enconder
        self.history = self.method_autoenconder.fit(X_train, X_train,
                                                    epochs=self.epochs,
                                                    batch_size=self.batch_size,
                                                    shuffle=True,
                                                    validation_data=(
                                                    X_validation, X_validation),
                                                    verbose=0)

    def transform(self, X):
        """ Template
        TO-DO: detailed explanation
        
        Parameters
        ----------


        

        """ 
        return self.method_enconder.predict(X)
    
def feature_learning(epochs, batch_size, name_dataset,
                     type_loss, value_encoding_dim,
                     X_train, X_test):
    """ Template
    TO-DO: detailed explanation

    Parameters
    ----------




    """ 
    
    autoEncoder_ = autoenconder(epochs=epochs,
                                batch_size=batch_size,
                                name_dataset=name_dataset,
                                type_loss=type_loss,
                                value_encoding_dim=value_encoding_dim)

    autoEncoder_.fit(X_train, X_test)

    X_train_encode = autoEncoder_.transform(X_train)
    X_test_encode = autoEncoder_.transform(X_test)

    return X_train_encode, X_test_encode, autoEncoder_