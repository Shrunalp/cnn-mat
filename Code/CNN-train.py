import wandb
wandb.login(key="<API Key>")
from wandb.keras import WandbMetricsLogger, WandbModelCheckpoint
import random
import numpy as np
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Conv2D, BatchNormalization, MaxPooling2D, Flatten, concatenate, Dropout
from tensorflow.keras.models import Model
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping
from Reformat_DREAM3D_Data import b2w

sz = "128"
#typ = "IPF_"
#typ = "FID_"
# ls = [typ+"BI_cell_data"+sz+".npy" , typ+"EI_BG_cell_data"+sz+".npy", typ+"EI_cell_data"+sz+".npy", 
#        typ+"EI_Cu_cell_data"+sz+".npy", typ+"RI_cell_data"+sz+".npy"]

typ = "PQ_"

ls = [typ+"BI_quat_data"+sz+".npy" , typ+"EI_BG_quat_data"+sz+".npy", typ+"EI_quat_data"+sz+".npy", 
        typ+"EI_Cu_quat_data"+sz+".npy", typ+"RI_quat_data"+sz+".npy"]
X = b2w(ls, 12)


sweep_config = {
    'method': "random",
    'metric': {
        'name': 'accuracy',
        'goal': 'maximize',
    },
    'parameters': {
        "optimizer": {
            "values": ['adam', 'sgd']
        },
        "conv1_filter": {
            "values": [8, 16, 32]
        },
        "conv2_filter": {
            "values": [8, 16, 32]
        },
        "layer_1": {
            "values": [128, 256, 512]
        },
        "layer_2": {
            "values": [64, 128, 256, 512]
        },
        "layer_3": {
            "values": [32, 64, 128, 256]
        },
        "dropout": {
            "values": [0.1, 0.2, 0.3, 0.5, 0.75]
        },
        "l1": {
            "values": [0.01, 0.1, 0.2, 0.25, 0.5, 0.75, 1]
        },
        "l2": {
            "values": [0.01, 0.1, 0.2, 0.25, 0.5, 0.75, 1]
        },
    },
}


sweep_id = wandb.sweep(sweep_config, project="EAGER_128")


def train():
    # Start a run, tracking hyperparameters
    wandb.init(
    # track hyperparameters and run metadata with wandb.config
        config={
            "conv1_filter": 8,
            "conv2_filter": 8,
            "dropout": random.uniform(0.01, 0.80),
            "layer_1": 128,
            "layer_2": 128, 
            "layer_3": 128, 
            "l1": 0.01, 
            "l2": 0.01, 
            "optimizer": "adam", 
            "metric": "accuracy",
        }
    )

    # Use wandb.config as your config
    config = wandb.config
    
    # Model Architecture
    input_shape=(64, 64, 4)

    nn_ph_model_input =  Input(shape=input_shape)

    x = Conv2D(config.conv1_filter, (3,3), strides=2, activation='relu', kernel_regularizer=tf.keras.regularizers.l1_l2(l1=config.l1, l2=config.l2), padding='same') (nn_ph_model_input)
    x = BatchNormalization() (x)
    x = MaxPooling2D((3,3), padding='same') (x)
    x = Dropout(config.dropout) (x)

    x = Conv2D(config.conv2_filter, (3,3), strides=2, activation='relu', kernel_regularizer=tf.keras.regularizers.l1_l2(l1=config.l1, l2=config.l2), padding='same') (x)
    x = BatchNormalization() (x)
    x = MaxPooling2D((3,3), padding='same') (x)
    x = Dropout(config.dropout) (x)


    x = Flatten()(x)

    x = Dense(config.layer_1, activation='relu') (x)
    x = Dense(config.layer_2, activation='relu') (x)
    x = Dense(config.layer_3, activation='relu') (x)
    x = Dropout(config.dropout) (x)

    output_layer = Dense(5, activation='softmax') (x)
    nn_ph_model = Model(nn_ph_model_input, output_layer)

    # Optimization
    nn_ph_model.compile(optimizer=config.optimizer, loss="categorical_crossentropy", metrics=[config.metric])
    
    # WandbMetricsLogger will log train and validation metrics to wandb
    # WandbModelCheckpoint will upload model checkpoints to wandb

    early_stopping = EarlyStopping(monitor='val_loss', patience=5, restore_best_weights=True)
    cb = [WandbMetricsLogger(log_freq=5), WandbModelCheckpoint("models"), early_stopping]

    history = nn_ph_model.fit(X[0], X[1], epochs= 50, validation_split=0.2, callbacks=cb)

wandb.agent(sweep_id, train, count=50)



