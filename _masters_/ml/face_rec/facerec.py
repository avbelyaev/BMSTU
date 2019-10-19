import os

import cv2
from PIL import Image
from numpy import asarray, savez_compressed, load
from mtcnn.mtcnn import MTCNN


FACES_FILE = 'faces.npz'


def extract_face(filename: str, required_size=(160, 160)):
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
    face_array = asarray(image)
    return face_array


def load_faces(directory):
    faces = []
    for filename in os.listdir(directory):
        face = extract_face(directory + filename)
        faces.append(face)
    return faces


def load_dataset(directory):
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
    return asarray(faces), asarray(labels)

def save_dataset_as_file(filename:str):
    train_faces, train_labels = load_dataset('data/train/')
    print(f'train faces/labels: {train_faces.shape}, {train_labels.shape}')

    test_faces, test_labels = load_dataset('data/test/')
    print(f'test faces/labels: {test_faces.shape}, {test_labels.shape}')

    # save whatever this fucks saves
    savez_compressed(filename, train_faces, train_labels, test_faces, test_labels)


def clusterize_from_file(filename: str):
    data = load(filename)
    train_faces, train_labels, test_faces, test_labels = data['arr_0'], data['arr_1'], data['arr_2'], data['arr_3']
    print('Loaded: ', train_faces.shape, train_labels.shape, test_faces.shape, test_labels.shape)


def main():
    file_with_extracted_faces_already_exists = os.path.isfile(FACES_FILE)
    if file_with_extracted_faces_already_exists:
        clusterize_from_file(FACES_FILE)
        return
    else:
        save_dataset_as_file(FACES_FILE)
    clusterize_from_file(FACES_FILE)


if __name__ == '__main__':
    main()
