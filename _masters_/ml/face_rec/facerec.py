import os

import cv2
from PIL import Image
import numpy as np
from mtcnn.mtcnn import MTCNN
from keras.models import load_model

FACES_FILE = 'faces.npz'
FACES_EMBEDDINGS_FILE = 'faces-embeddings.npz'


# https://machinelearningmastery.com/how-to-develop-a-face-recognition-system-using-facenet-in-keras-and-an-svm-classifier/


def extract_face(filename: str, required_size=(160, 160)) -> np.ndarray:
    img = cv2.cvtColor(cv2.imread(filename), cv2.COLOR_BGR2RGB)

    detector = MTCNN()
    results = detector.detect_faces(img)
    # extract the bounding box from the first face
    x1, y1, width, height = results[0]['box']
    x1, y1 = abs(x1), abs(y1)
    x2, y2 = x1 + width, y1 + height

    face = img[y1:y2, x1:x2]
    # resize pixels to the model size
    # TODO replace this shit with opencv
    image = Image.fromarray(face)
    image = image.resize(required_size)
    face_array = np.asarray(image)
    return face_array


def load_faces(directory: str) -> list:
    faces = []
    for filename in os.listdir(directory):
        face = extract_face(directory + filename)
        faces.append(face)
    return faces


def load_dataset(directory: str) -> (list, list):
    faces, labels = [], []
    for subdir in os.listdir(directory):
        path = directory + subdir + '/'
        if not os.path.isdir(path):
            continue
        # load all faces in the subdirectory
        faces_loaded = load_faces(path)
        labels_loaded = [subdir for _ in range(len(faces_loaded))]
        print('>loaded %d examples for class: %s' % (len(faces_loaded), subdir))

        faces.extend(faces_loaded)
        labels.extend(labels_loaded)
    return np.asarray(faces), np.asarray(labels)


def save_dataset_as_file(filename: str):
    train_faces, train_labels = load_dataset('data/train/')
    print(f'train faces/labels: {train_faces.shape}, {train_labels.shape}')

    test_faces, test_labels = load_dataset('data/test/')
    print(f'test faces/labels: {test_faces.shape}, {test_labels.shape}')

    # save whatever this fucks saves
    np.savez_compressed(filename, train_faces, train_labels, test_faces, test_labels)


# get the face embedding for one face
def get_embedding(model, face_pixels):
    # scale pixel values
    face_pixels = face_pixels.astype('float32')
    # standardize pixel values across channels (global)
    mean, std = face_pixels.mean(), face_pixels.std()
    face_pixels = (face_pixels - mean) / std
    # transform face into one sample
    samples = np.expand_dims(face_pixels, axis=0)
    # make prediction to get embedding
    yhat = model.predict(samples)
    return yhat[0]


def save_embeddings_as_file():
    data = np.load(FACES_FILE)
    train_faces, train_labels, test_faces, test_labels = data['arr_0'], data['arr_1'], data['arr_2'], data['arr_3']
    print('Loaded: ', train_faces.shape, train_labels.shape, test_faces.shape, test_labels.shape)

    # load MOAR of that shit!
    model = load_model('facenet_keras.h5')

    newTrainX = list()
    for face_pixels in train_faces:
        embedding = get_embedding(model, face_pixels)
        newTrainX.append(embedding)
    newTrainX = np.asarray(newTrainX)
    print(newTrainX.shape)
    # convert each face in the test set to an embedding
    newTestX = list()
    for face_pixels in test_faces:
        embedding = get_embedding(model, face_pixels)
        newTestX.append(embedding)
    newTestX = np.asarray(newTestX)
    print(newTestX.shape)

    np.savez_compressed(FACES_EMBEDDINGS_FILE, newTrainX, train_labels, newTestX, test_labels)


def main():
    file_with_extracted_faces_already_exists = os.path.isfile(FACES_FILE)
    if file_with_extracted_faces_already_exists:
        save_embeddings_as_file()

        return
    else:
        save_dataset_as_file(FACES_FILE)

    save_embeddings_as_file()


if __name__ == '__main__':
    main()
