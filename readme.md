# Numerical Mathematics II for Engineers – Homework Repository (WiSe 2025/26)

**Group 9**

This repository contains the homework solutions for **Numerical Mathematics II for Engineers** in the Winter Semester 2025/26.

## Repo Structure

- `homework#/` – Each contains homework solutions of homework #

## How to convert .ipynb to .pdf

```
# install jupyter and pandoc
brew install pandoc
brew install --cask mactex

curl -sL "https://yihui.org/tinytex/install-bin-unix.sh" | sh

sudo /usr/local/bin/tlmgr path add
xelatex --version
tlmgr install collection-latexextra collection-fontsrecommended collection-latexrecommended
```

```
jupyter nbconvert --to pdf [notebook].ipynb
```

## How to Submit the homework

Hand in the solution in one folder labeled hw1_group[group number] and containing:

- One pdf for the theoretical questions and comments on the numerical results
- One python file per programming exercise.

Write the group number and all names of your members in each file.

## How to Run

1. Clone the repository:

```bash
git clone https://github.com/TU-Berlin-SC/Numerical-Mathematics-II-for-Engineers.git
cd Numerical-Mathematics-II-for-Engineers
```

## Homework Lists

|     | hw         | solution                           | Date       |
| --- | ---------- | ---------------------------------- | ---------- |
| 0   | Homework 0 | [Homework0 solution](./homework0/) | 20.10.2025 |
| 1   | Homework 1 | [Homework1 solution](./homework1/) | 22.10.2025 |
|     |            |                                    |            |
