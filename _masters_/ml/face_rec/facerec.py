from os import listdir

import cv2
from PIL import Image
from matplotlib import pyplot
from numpy import asarray
from mtcnn.mtcnn import MTCNN


def extract_face(filename: str, required_size=(160, 160)):
    img = cv2.cvtColor(cv2.imread(filename), cv2.COLOR_BGR2RGB)
    # image = Image.open(filename)
    # image = image.convert('RGB')
    # pixels = asarray(image)

    detector = MTCNN()
    results = detector.detect_faces(img)
    # extract the bounding box from the first face
    x1, y1, width, height = results[0]['box']
    # bug fix
    x1, y1 = abs(x1), abs(y1)
    x2, y2 = x1 + width, y1 + height
    # extract the face
    face = img[y1:y2, x1:x2]
    # resize pixels to the model size
    image = Image.fromarray(face)
    image = image.resize(required_size)
    face_array = asarray(image)
    return face_array


def main():
    folder = 'data/LeBron_James/'
    i = 1
    for filename in listdir(folder):
        face = extract_face(folder + filename)
        print(i, face.shape)
        # plot
        pyplot.subplot(2, 7, i)
        pyplot.axis('off')
        pyplot.imshow(face)
        i += 1
    pyplot.show()


if __name__ == '__main__':
    main()
