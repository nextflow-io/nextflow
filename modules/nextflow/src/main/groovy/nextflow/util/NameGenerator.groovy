/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.util

import groovy.transform.CompileStatic
import groovy.transform.PackageScope

/**
 * Generate a random mnemonic name
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class NameGenerator {

    @PackageScope
    static List<String> ADJECTIVES = [

            "admiring",
            "adoring",
            "agitated",
            "amazing",
            "angry",
            "astonishing",
            "awesome",
            "backstabbing",
            "berserk",
            "big",
            "boring",
            "chaotic",
            "cheeky",
            "cheesy",
            "clever",
            "compassionate",
            "condescending",
            "confident",
            "cranky",
            "crazy",
            "curious",
            "deadly",
            "desperate",
            "determined",
            "distracted",
            "distraught",
            "disturbed",
            "dreamy",
            "drunk",
            "ecstatic",
            "elated",
            "elegant",
            "evil",
            "exotic",
            "extravagant",
            "fabulous",
            "fervent",
            "festering",
            "focused",
            "friendly",
            "furious",
            "gigantic",
            "gloomy",
            "golden",
            "goofy",
            "grave",
            "happy",
            "high",
            "hopeful",
            "hungry",
            "infallible",
            "insane",
            "intergalactic",
            "irreverent",
            "jolly",
            "jovial",
            "kickass",
            "lethal",
            "lonely",
            "loving",
            "mad",
            "magical",
            "maniac",
            "marvelous",
            "mighty",
            "modest",
            "nasty",
            "naughty",
            "nauseous",
            "nice",
            "nostalgic",
            "peaceful",
            "pedantic",
            "pensive",
            "prickly",
            "reverent",
            "ridiculous",
            "romantic",
            "sad",
            "scruffy",
            "serene",
            "sharp",
            "shrivelled",
            "sick",
            "silly",
            "sleepy",
            "small",
            "soggy",
            "special",
            "spontaneous",
            "stoic",
            "stupefied",
            "suspicious",
            "tender",
            "thirsty",
            "tiny",
            "trusting",
            "voluminous",
            "wise",
            "zen"

    ]

    @PackageScope
    static final List<String> NAMES = [

            // Maria Gaetana Agnesi - Italian mathematician, philosopher, theologian and humanitarian. She was the first woman to write a mathematics handbook and the first woman appointed as a Mathematics Professor at a University. https://en.wikipedia.org/wiki/Maria_Gaetana_Agnesi
            "agnesi",

            // Muhammad ibn Jābir al-Ḥarrānī al-Battānī was a founding father of astronomy. https://en.wikipedia.org/wiki/Al-Battani
            "albattani",

            // Frances E. Allen, became the first female IBM Fellow in 1989. In 2006, she became the first female recipient of the ACM's Turing Award. https://en.wikipedia.org/wiki/Frances_E._Allen
            "allen",

            // June Almeida - Scottish virologist who took the first pictures of the rubella virus. https://en.wikipedia.org/wiki/June_Almeida
            "almeida",

            // André-Marie Ampère - French physicist and mathematician, one of the founders of the science of classical electromagnetism. https://en.wikipedia.org/wiki/Andr%C3%A9-Marie_Amp%C3%A8re
            "ampere",

            // Archimedes was a physicist, engineer and mathematician who invented too many things to list them here. https://en.wikipedia.org/wiki/Archimedes
            "archimedes",

            // Maria Ardinghelli - Italian translator, mathematician and physicist. https://en.wikipedia.org/wiki/Maria_Ardinghelli
            "ardinghelli",

            // Aryabhata - Ancient Indian mathematician-astronomer during 476-550 CE. https://en.wikipedia.org/wiki/Aryabhata
            "aryabhata",

            // Wanda Austin - Wanda Austin is the President and CEO of The Aerospace Corporation, a leading architect for the US security space programs. https://en.wikipedia.org/wiki/Wanda_Austin
            "austin",

            // Lorenzo Romano Amedeo Carlo Avogadro was an Italian scientist, most noted for his contribution to molecular theory now known as Avogadro's law. https://en.wikipedia.org/wiki/Amedeo_Avogadro
            "avogadro",

            // Charles Babbage invented the concept of a programmable computer. https://en.wikipedia.org/wiki/Charles_Babbage.
            "babbage",

            // Leo Baekeland - Belgian-American scientist known for inventing Bakelite. Considered the "Father of the Plastics Industry". https://en.wikipedia.org/wiki/Leo_Baekeland
            "baekeland",

            // Stefan Banach - Polish mathematician, was one of the founders of modern functional analysis. https://en.wikipedia.org/wiki/Stefan_Banach
            "banach",

            // John Bardeen co-invented the transistor. https://en.wikipedia.org/wiki/John_Bardeen
            "bardeen",

            // Jean Jennings Bartik, born Betty Jean Jennings - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Jean_Bartik
            "bartik",

            // Laura Bassi, the world's first female professor. https://en.wikipedia.org/wiki/Laura_Bassi
            "bassi",

            // Antoine César Becquerel, Alexandre-Edmond Becquerel, Antoine Henri Becquerel and Jean Becquerel... Lineage of French scientists, working on electric and luminescent phenomena, solar spectrum, magnetism, electricity, optics, radioactivity... https://en.wikipedia.org/wiki/Becquerel_(disambiguation)
            "becquerel",

            // Alexander Graham Bell - an eminent Scottish-born scientist, inventor, engineer and innovator who is credited with inventing the first practical telephone. https://en.wikipedia.org/wiki/Alexander_Graham_Bell
            "bell",

            // Claude Bernard - French physiologist, was one of the first to suggest the use of blind experiments to ensure the objectivity of scientific observations. https://en.wikipedia.org/wiki/Claude_Bernard
            "bernard",

            // Homi J Bhabha - was an Indian nuclear physicist, founding director, and professor of physics at the Tata Institute of Fundamental Research. Colloquially known as "father of Indian nuclear programme". https://en.wikipedia.org/wiki/Homi_J._Bhabha
            "bhabha",

            // Bhaskara II - Ancient Indian mathematician-astronomer whose work on calculus predates Newton and Leibniz by over half a millennium. https://en.wikipedia.org/wiki/Bh%C4%81skara_II#Calculus
            "bhaskara",

            // Elizabeth Blackwell - American doctor and first American woman to receive a medical degree. https://en.wikipedia.org/wiki/Elizabeth_Blackwell
            "blackwell",

            // Niels Bohr is the father of quantum theory. https://en.wikipedia.org/wiki/Niels_Bohr.
            "bohr",

            // Kathleen Booth, she's credited with writing the first assembly language. https://en.wikipedia.org/wiki/Kathleen_Booth
            "booth",

            // Anita Borg - Anita Borg was the founding director of the Institute for Women and Technology (IWT). https://en.wikipedia.org/wiki/Anita_Borg
            "borg",

            // Satyendra Nath Bose - He provided the foundation for Bose–Einstein statistics and the theory of the Bose–Einstein condensate. https://en.wikipedia.org/wiki/Satyendra_Nath_Bose
            "bose",

            // Evelyn Boyd Granville - She was one of the first African-American woman to receive a Ph.D. in mathematics; she earned it in 1949 from Yale University. https://en.wikipedia.org/wiki/Evelyn_Boyd_Granville
            "boyd",

            // Brahmagupta - Ancient Indian mathematician during 598-670 CE who gave rules to compute with zero. https://en.wikipedia.org/wiki/Brahmagupta#Zero
            "brahmagupta",

            // Walter Houser Brattain co-invented the transistor. https://en.wikipedia.org/wiki/Walter_Houser_Brattain
            "brattain",

            // Sydney Brenner - South African biologist who worked on the genetic code, and other areas of molecular biology. https://en.wikipedia.org/wiki/Sydney_Brenner
            "brenner",

            // Emmett Brown invented time travel. https://en.wikipedia.org/wiki/Emmett_Brown (thanks Brian Goff)
            "brown",

            // Santiago Ramón y Cajal - Spanish, pathologist, histologist, neuroscientist, is considered the father of neuroscience. https://en.wikipedia.org/wiki/Santiago_Ram%C3%B3n_y_Cajal
            "cajal",

            // Georg Ferdinand Ludwig Philipp Cantor - German mathematician, invented set theory, which has become a fundamental theory in mathematics. https://en.wikipedia.org/wiki/Georg_Cantor
            "cantor",

            // Michelangelo Merisi da Caravaggio was an Italian painter active. https://en.wikipedia.org/wiki/Caravaggio
            "caravaggio",

            // Albertina Carlsson - Swedish zoologist, referred to as the first Swedish woman to have performed scientific studies in zoology. https://en.wikipedia.org/wiki/Albertina_Carlsson
            "carlsson",

            // Rachel Carson - American marine biologist and conservationist, her book Silent Spring and other writings are credited with advancing the global environmental movement. https://en.wikipedia.org/wiki/Rachel_Carson
            "carson",

            // Anders Celsius - Swedish astronomer, physicist and mathematician, proposed the Celsius temperature scale which bears his name. https://en.wikipedia.org/wiki/Anders_Celsius
            "celsius",

            // Subrahmanyan Chandrasekhar - Astrophysicist known for his mathematical theory on different stages and evolution in structures of the stars. He has won nobel prize for physics. https://en.wikipedia.org/wiki/Subrahmanyan_Chandrasekhar
            "chandrasekhar",

            // George M. Church - American geneticist, molecular engineer, and chemist. https://en.wikipedia.org/wiki/George_M._Church
            "church",

            // Jane Colden - American botanist widely considered the first female American botanist. https://en.wikipedia.org/wiki/Jane_Colden
            "colden",

            // Gerty Theresa Cori - American biochemist who became the third woman—and first American woman—to win a Nobel Prize in science, and the first woman to be awarded the Nobel Prize in Physiology or Medicine. Cori was born in Prague. https://en.wikipedia.org/wiki/Gerty_Cori
            "cori",

            // Charles-Augustin de Coulomb - French physicist, developed Coulomb's law, the definition of the electrostatic force of attraction and repulsion, but also did important work on friction. https://en.wikipedia.org/wiki/Charles-Augustin_de_Coulomb
            "coulomb",

            // Seymour Roger Cray was an American electrical engineer and supercomputer architect who designed a series of computers that were the fastest in the world for decades. https://en.wikipedia.org/wiki/Seymour_Cray
            "cray",

            // Francis Crick - British molecular biologist, biophysicist, and neuroscientist, most noted for being a co-discoverer of the structure of the DNA molecule. https://en.wikipedia.org/wiki/Francis_Crick
            "crick",

            // This entry reflects a husband and wife team who worked together:
            // Marie Curie discovered radioactivity. https://en.wikipedia.org/wiki/Marie_Curie
            // Pierre Curie was a pioneer in crystallography, magnetism, piezoelectricity and radioactivity. https://en.wikipedia.org/wiki/Pierre_Curie
            "curie",

            // This entry reflects a husband and wife team who worked together:
            // Joan Curran was a Welsh scientist who developed radar and invented chaff, a radar countermeasure. https://en.wikipedia.org/wiki/Joan_Curran
            // Samuel Curran was an Irish physicist who worked alongside his wife during WWII and invented the proximity fuse. https://en.wikipedia.org/wiki/Samuel_Curran
            "curran",

            // Georges Cuvier - French naturalist and zoologist, was instrumental in establishing the fields of comparative anatomy and paleontology through his work in comparing living animals with fossils. https://en.wikipedia.org/wiki/Georges_Cuvier
            "cuvier",

            // Jean-Baptiste le Rond d'Alembert - French mathematician, mechanician, physicist, philosopher, and music theorist. https://en.wikipedia.org/wiki/Jean_le_Rond_d%27Alembert
            "dalembert",

            // Charles Darwin established the principles of natural evolution. https://en.wikipedia.org/wiki/Charles_Darwin.
            "darwin",

            // Leonardo Da Vinci invented too many things to list here. https://en.wikipedia.org/wiki/Leonardo_da_Vinci.
            "davinci",

            // René Descartes - French philosopher, mathematician, and scientist. His best known philosophical statement is "Cogito ergo sum". https://en.wikipedia.org/wiki/Ren%C3%A9_Descartes
            "descartes",

            // Edsger Wybe Dijkstra was a Dutch computer scientist and mathematical scientist. https://en.wikipedia.org/wiki/Edsger_W._Dijkstra.
            "dijkstra",

            // Donna Dubinsky - played an integral role in the development of personal digital assistants (PDAs) serving as CEO of Palm, Inc. and co-founding Handspring. https://en.wikipedia.org/wiki/Donna_Dubinsky
            "dubinsky",

            // Annie Easley - She was a leading member of the team which developed software for the Centaur rocket stage and one of the first African-Americans in her field. https://en.wikipedia.org/wiki/Annie_Easley
            "easley",

            // Thomas Alva Edison, prolific inventor. https://en.wikipedia.org/wiki/Thomas_Edison
            "edison",

            // Albert Einstein invented the general theory of relativity. https://en.wikipedia.org/wiki/Albert_Einstein
            "einstein",

            // Eva Ekeblad - Swedish countess who was a salon hostess, agronomist, and scientist, discovered a method in 1746 to make alcohol and flour from potatoes. She was the first female member of the Royal Swedish Academy of Sciences. https://en.wikipedia.org/wiki/Eva_Ekeblad
            "ekeblad",

            // Gertrude Elion - American biochemist, pharmacologist and the 1988 recipient of the Nobel Prize in Medicine. https://en.wikipedia.org/wiki/Gertrude_Elion
            "elion",

            // Douglas Engelbart gave the mother of all demos: https://en.wikipedia.org/wiki/Douglas_Engelbart
            "engelbart",

            // Euclid invented geometry. https://en.wikipedia.org/wiki/Euclid
            "euclid",

            // Leonhard Euler invented large parts of modern mathematics. https://de.wikipedia.org/wiki/Leonhard_Euler
            "euler",

            // Federico Faggin is an Italian physicist, inventor and entrepreneur, widely known for designing the first commercial microprocessor. https://en.wikipedia.org/wiki/Federico_Faggin
            "faggin",

            // Pierre de Fermat pioneered several aspects of modern mathematics. https://en.wikipedia.org/wiki/Pierre_de_Fermat
            "fermat",

            // Enrico Fermi invented the first nuclear reactor. https://en.wikipedia.org/wiki/Enrico_Fermi.
            "fermi",

            // Richard Feynman was a key contributor to quantum mechanics and particle physics. https://en.wikipedia.org/wiki/Richard_Feynman
            "feynman",

            // Jean-Baptiste Joseph Fourier - a French mathematician and physicist, best known for initiating the investigation of Fourier series and their applications to problems of heat transfer and vibrations, is also generally credited with the discovery of the greenhouse effect. https://en.wikipedia.org/wiki/Joseph_Fourier
            "fourier",

            // This entry reflects two completely unrelated scientists:
            // Benjamin Franklin is famous for his experiments in electricity and the invention of the lightning rod. https://en.wikipedia.org/wiki/Benjamin_Franklin
            // Rosalind Franklin - English chemist and X-ray crystallographer whose research was critical to the understanding of DNA. https://en.wikipedia.org/wiki/Rosalind_Franklin
            "franklin",

            // Galileo was a founding father of modern astronomy, and faced politics and obscurantism to establish scientific truth.  https://en.wikipedia.org/wiki/Galileo_Galilei
            "galileo",

            // William Henry "Bill" Gates III is an American business magnate, philanthropist, investor, computer programmer, and inventor. https://en.wikipedia.org/wiki/Bill_Gates
            "gates",

            // Marthe Gautier discovered the link of diseases to chromosome abnormalities. https://en.wikipedia.org/wiki/Marthe_Gautier
            "gautier",

            // Walter Gilbert is an American biochemist, developed a DNA sequencing method and first proposed the existence of introns and exons. https://en.wikipedia.org/wiki/Walter_Gilbert
            "gilbert",

            // Adele Goldberg, was one of the designers and developers of the Smalltalk language. https://en.wikipedia.org/wiki/Adele_Goldberg_(computer_scientist)
            "goldberg",

            // Adele Goldstine, born Adele Katz, wrote the complete technical description for the first electronic digital computer, ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Adele_Goldstine
            "goldstine",

            // Shafi Goldwasser is a computer scientist known for creating theoretical foundations of modern cryptography. Winner of 2012 ACM Turing Award. https://en.wikipedia.org/wiki/Shafi_Goldwasser
            "goldwasser",

            // James Golick, all around gangster.
            "golick",

            // Jane Goodall - British primatologist, ethologist, and anthropologist who is considered to be the world's foremost expert on chimpanzees. https://en.wikipedia.org/wiki/Jane_Goodall
            "goodall",

            // Johannes Gensfleisch zur Laden zum Gutenberg - Inventor of movable-type printing press. https://en.wikipedia.org/wiki/Johannes_Gutenberg
            "gutenberg",

            // Margaret Hamilton - Director of the Software Engineering Division of the MIT Instrumentation Laboratory, which developed on-board flight software for the Apollo space program. https://en.wikipedia.org/wiki/Margaret_Hamilton_(scientist)
            "hamilton",

            // Stephen Hawking pioneered the field of cosmology by combining general relativity and quantum mechanics. https://en.wikipedia.org/wiki/Stephen_Hawking
            "hawking",

            // Werner Heisenberg was a founding father of quantum mechanics. https://en.wikipedia.org/wiki/Werner_Heisenberg
            "heisenberg",

            // Jaroslav Heyrovský was the inventor of the polarographic method, father of the electroanalytical method, and recipient of the Nobel Prize in 1959. His main field of work was polarography. https://en.wikipedia.org/wiki/Jaroslav_Heyrovsk%C3%BD
            "heyrovsky",

            // David Hilbert - German mathematician, recognized as one of the most influential and universal mathematicians of the 19th and early 20th centuries. he is known as one of the founders of proof theory and mathematical logic, as well as for being among the first to distinguish between mathematics and metamathematics. https://en.wikipedia.org/wiki/David_Hilbert
            "hilbert",

            // Dorothy Hodgkin was a British biochemist, credited with the development of protein crystallography. She was awarded the Nobel Prize in Chemistry in 1964. https://en.wikipedia.org/wiki/Dorothy_Hodgkin
            "hodgkin",

            // Erna Schneider Hoover revolutionized modern communication by inventing a computerized telephone switching method. https://en.wikipedia.org/wiki/Erna_Schneider_Hoover
            "hoover",

            // Grace Hopper developed the first compiler for a computer programming language and  is credited with popularizing the term "debugging" for fixing computer glitches. https://en.wikipedia.org/wiki/Grace_Hopper
            "hopper",

            // Frances Hugle, she was an American scientist, engineer, and inventor who contributed to the understanding of semiconductors, integrated circuitry, and the unique electrical principles of microscopic materials. https://en.wikipedia.org/wiki/Frances_Hugle
            "hugle",

            // Hypatia - Greek Alexandrine Neoplatonist philosopher in Egypt who was one of the earliest mothers of mathematics. https://en.wikipedia.org/wiki/Hypatia
            "hypatia",

            // Yeong-Sil Jang was a Korean scientist and astronomer during the Joseon Dynasty; he invented the first metal printing press and water gauge. https://en.wikipedia.org/wiki/Jang_Yeong-sil
            "jang",

            // Jean Jennings Bartik, born Betty Jean Jennings - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Jean_Bartik
            "jennings",

            // Mary Lou Jepsen, was the founder and chief technology officer of One Laptop Per Child (OLPC), and the founder of Pixel Qi. https://en.wikipedia.org/wiki/Mary_Lou_Jepsen
            "jepsen",

            // Irène Joliot-Curie - French scientist who was awarded the Nobel Prize for Chemistry in 1935. Daughter of Marie and Pierre Curie. https://en.wikipedia.org/wiki/Ir%C3%A8ne_Joliot-Curie
            "joliot",

            // Karen Spärck Jones came up with the concept of inverse document frequency, which is used in most search engines today. https://en.wikipedia.org/wiki/Karen_Sp%C3%A4rck_Jones
            "jones",

            // A. P. J. Abdul Kalam - is an Indian scientist aka Missile Man of India for his work on the development of ballistic missile and launch vehicle technology. https://en.wikipedia.org/wiki/A._P._J._Abdul_Kalam
            "kalam",

            // Rudolf E. Kálmán, American engineer and mathematican of Hungarian origin. One of the inventors of the smoother/predictor commonly known as "Kalman Filter". https://en.wikipedia.org/wiki/Rudolf_E._K%C3%A1lm%C3%A1n
            "kalman",

            // Susan Kare, created the icons and many of the interface elements for the original Apple Macintosh in the 1980s, and was an original employee of NeXT, working as the Creative Director. https://en.wikipedia.org/wiki/Susan_Kare
            "kare",

            // Mary Kenneth Keller, Sister Mary Kenneth Keller became the first American woman to earn a PhD in Computer Science in 1965. https://en.wikipedia.org/wiki/Mary_Kenneth_Keller
            "keller",

            // Har Gobind Khorana - Indian-American biochemist who shared the 1968 Nobel Prize for Physiology. https://en.wikipedia.org/wiki/Har_Gobind_Khorana
            "khorana",

            // Jack Kilby invented silicone integrated circuits and gave Silicon Valley its name. https://en.wikipedia.org/wiki/Jack_Kilby
            "kilby",

            // Motoo Kimura - Japanese biologist introduced the neutral theory of molecular evolution. https://en.wikipedia.org/wiki/Motoo_Kimura
            "kimura",

            // Maria Kirch - German astronomer and first woman to discover a comet. https://en.wikipedia.org/wiki/Maria_Margarethe_Kirch
            "kirch",

            // Donald Knuth - American computer scientist, author of "The Art of Computer Programming" and creator of the TeX typesetting system. https://en.wikipedia.org/wiki/Donald_Knuth
            "knuth",

            // Robert Kock - German physician, considered the founder of modern bacteriology. https://en.wikipedia.org/wiki/Robert_Koch
            "koch",

            // Sophie Kowalevski - Russian mathematician responsible for important original contributions to analysis, differential equations and mechanics. https://en.wikipedia.org/wiki/Sofia_Kovalevskaya
            "kowalevski",

            // Marie-Jeanne de Lalande - French astronomer, mathematician and cataloguer of stars. https://en.wikipedia.org/wiki/Marie-Jeanne_de_Lalande
            "lalande",

            // Jean-Baptiste Lamarck - French naturalist. He was a soldier, biologist, academic, and an early proponent of the idea that biological evolution occurred and proceeded in accordance with natural laws. https://en.wikipedia.org/wiki/Jean-Baptiste_Lamarck
            "lamarck",

            // Hedy Lamarr - Actress and inventor. The principles of her work are now incorporated into modern Wi-Fi, CDMA and Bluetooth technology. https://en.wikipedia.org/wiki/Hedy_Lamarr
            "lamarr",

            // Leslie B. Lamport - American computer scientist. Lamport is best known for his seminal work in distributed systems and was the winner of the 2013 Turing Award. https://en.wikipedia.org/wiki/Leslie_Lamport
            "lamport",

            // Pierre-Simon, marquis de Laplace, French scholar whose work was important to the development of mathematics, statistics, physics and astronomy. https://en.wikipedia.org/wiki/Pierre-Simon_Laplace
            "laplace",

            // Antoine Lavoisier - French chemist central to the 18th-century chemical revolution and had a large influence on both the history of chemistry and the history of biology. https://en.wikipedia.org/wiki/Antoine_Lavoisier
            "lavoisier",

            // Mary Leakey - British paleoanthropologist who discovered the first fossilized Proconsul skull. https://en.wikipedia.org/wiki/Mary_Leakey
            "leakey",

            // Henrietta Swan Leavitt - she was an American astronomer who discovered the relation between the luminosity and the period of Cepheid variable stars. https://en.wikipedia.org/wiki/Henrietta_Swan_Leavitt
            "leavitt",

            // Guillaume Joseph Hyacinthe Jean-Baptiste Le Gentil de la Galaisière - French astronomer, part of the international collaborative project to measure the distance to the Sun, by observing the transit of Venus at different points on the earth. Mainly known for being one of the most unfortunate and unlucky scientist ever. https://en.wikipedia.org/wiki/Guillaume_Le_Gentil
            "legentil",

            // Ruth Lichterman - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Ruth_Teitelbaum
            "lichterman",

            // Carl Linnaeus (also known as Carl von Linné, Carolus Linnæus or Carolus a Linné) - Swedish botanist, physician, and zoologist, father of modern taxonomy, he formalised the modern system of naming organisms called binomial nomenclature. https://en.wikipedia.org/wiki/Carl_Linnaeus
            "linnaeus",

            // Barbara Liskov - co-developed the Liskov substitution principle. Liskov was also the winner of the Turing Prize in 2008. https://en.wikipedia.org/wiki/Barbara_Liskov
            "liskov",

            // Konrad Zacharias Lorenz - was an Austrian zoologist, ethologist, and ornithologist founder of the modern ethology, study of animal behavior. https://en.wikipedia.org/wiki/Konrad_Lorenz
            "lorenz",

            // Ada Lovelace invented the first algorithm. https://en.wikipedia.org/wiki/Ada_Lovelace (thanks James Turnbull)
            "lovelace",

            // Auguste and Louis Lumière - the first filmmakers in history. https://en.wikipedia.org/wiki/Auguste_and_Louis_Lumi%C3%A8re
            "lumiere",

            // René François Ghislain Magritte - Belgian surrealist artist. https://en.wikipedia.org/wiki/René_Magritte
            "magritte",

            // Mahavira - Ancient Indian mathematician during 9th century AD who discovered basic algebraic identities. https://en.wikipedia.org/wiki/Mah%C4%81v%C4%ABra_(mathematician)
            "mahavira",

            // Ettore Majorana was an Italian theoretical physicist who worked on neutrino masses. https://en.wikipedia.org/wiki/Ettore_Majorana
            "majorana",

            // Benoit B. Mandelbrot - Polish-born, French and American mathematician, recognized for his contribution to the field of fractal geometry, which included coining the word "fractal". https://en.wikipedia.org/wiki/Benoit_Mandelbrot
            "mandelbrot",

            // Guglielmo Marconi, 1st Marquis of Marconi was an Italian inventor and electrical engineer known for his pioneering work on long-distance radio transmission. https://en.wikipedia.org/wiki/Guglielmo_Marconi
            "marconi",

            // Maria Mayer - American theoretical physicist and Nobel laureate in Physics for proposing the nuclear shell model of the atomic nucleus. https://en.wikipedia.org/wiki/Maria_Mayer
            "mayer",

            // John McCarthy invented LISP. https://en.wikipedia.org/wiki/John_McCarthy_(computer_scientist)
            "mccarthy",

            // Barbara McClintock - a distinguished American cytogeneticist, 1983 Nobel Laureate in Physiology or Medicine for discovering transposons. https://en.wikipedia.org/wiki/Barbara_McClintock
            "mcclintock",

            // Malcolm McLean invented the modern shipping container: https://en.wikipedia.org/wiki/Malcom_McLean
            "mclean",

            // Kay McNulty - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Kathleen_Antonelli
            "mcnulty",

            // Lise Meitner - Austrian/Swedish physicist who was involved in the discovery of nuclear fission. The element meitnerium is named after her. https://en.wikipedia.org/wiki/Lise_Meitner
            "meitner",

            // Gregor Mendel - Scientist and Augustinian friar. Gained posthumous recognition as the founder of the modern science of genetics. https://en.wikipedia.org/wiki/Gregor_Mendel
            "mendel",

            // Carla Meninsky, was the game designer and programmer for Atari 2600 games Dodge 'Em and Warlords. https://en.wikipedia.org/wiki/Carla_Meninsky
            "meninsky",

            // Gerardus Mercator - German-Flemish cartographer, geographer and cosmographer. Known for creating a 1569 world map based on rhumb lines. https://en.wikipedia.org/wiki/Gerardus_Mercator
            "mercator",

            // Johanna Mestorf - German prehistoric archaeologist and first female museum director in Germany. https://en.wikipedia.org/wiki/Johanna_Mestorf
            "mestorf",

            // Antonio Santi Giuseppe Meucci was an Italian inventor best known for developing a voice-communication apparatus that several sources credit as the first telephone. https://en.wikipedia.org/wiki/Antonio_Meucci
            "meucci",

            // Friedrich Miescher - Swiss physician and biologist. He was the first researcher to isolate nucleic acid. https://en.wikipedia.org/wiki/Friedrich_Miescher
            "miescher",

            // Marvin Minsky - Pioneer in Artificial Intelligence, co-founder of the MIT's AI Lab, won the Turing Award in 1969. https://en.wikipedia.org/wiki/Marvin_Minsky
            "minsky",

            // Maryam Mirzakhani - an Iranian mathematician and the first woman to win the Fields Medal. https://en.wikipedia.org/wiki/Maryam_Mirzakhani
            "mirzakhani",

            // Jacques Lucien Monod - A french molecular biologist who studied genetic control of enzyme and virus synthesis using as model the lac operon of E. coli. https://en.wikipedia.org/wiki/Jacques_Monod
            "monod",

            // Rita Levi-Montalcini - Won Nobel Prize in Physiology or Medicine jointly with colleague Stanley Cohen for the discovery of nerve growth factor (https://en.wikipedia.org/wiki/Rita_Levi-Montalcini)
            "montalcini",

            // Samuel Morse - contributed to the invention of a single-wire telegraph system based on European telegraphs and was a co-developer of the Morse code. https://en.wikipedia.org/wiki/Samuel_Morse
            "morse",

            // Ian Murdock - founder of the Debian project. https://en.wikipedia.org/wiki/Ian_Murdock
            "murdock",

            // Isaac Newton invented classic mechanics and modern optics. https://en.wikipedia.org/wiki/Isaac_Newton
            "newton",

            // Florence Nightingale, more prominently known as a nurse, was also the first female member of the Royal Statistical Society and a pioneer in statistical graphics. https://en.wikipedia.org/wiki/Florence_Nightingale#Statistics_and_sanitary_reform
            "nightingale",

            // Alfred Nobel - a Swedish chemist, engineer, innovator, and armaments manufacturer (inventor of dynamite). https://en.wikipedia.org/wiki/Alfred_Nobel
            "nobel",

            // Emmy Noether, German mathematician. Noether's Theorem is named after her. https://en.wikipedia.org/wiki/Emmy_Noether
            "noether",

            // Frances 'Poppy' Northcutt. Poppy Northcutt was the first woman to work as part of NASA’s Mission Control. https://en.wikipedia.org/wiki/Frances_Northcutt - http://www.businessinsider.com/poppy-northcutt-helped-apollo-astronauts-2014-12
            "northcutt",

            // Robert Noyce invented silicone integrated circuits and gave Silicon Valley its name. https://en.wikipedia.org/wiki/Robert_Noyce
            "noyce",

            // Pāṇini - Ancient Indian linguist and grammarian from 4th century CE who worked on the world's first formal system. https://en.wikipedia.org/wiki/P%C4%81%E1%B9%87ini#Comparison_with_modern_formal_systems
            "panini",

            // Ambroise Paré invented modern surgery. https://en.wikipedia.org/wiki/Ambroise_Par%C3%A9
            "pare",

            // Louis Pasteur discovered vaccination, fermentation and pasteurization. https://en.wikipedia.org/wiki/Louis_Pasteur.
            "pasteur",

            // Linus Carl Pauling - American biochemist considered one of the founders of the fields of quantum chemistry and molecular biology. He was also a peace activist and was awarded the Nobel Peace Prize. https://en.wikipedia.org/wiki/Linus_Pauling
            "pauling",

            // Cecilia Payne-Gaposchkin was an astronomer and astrophysicist who, in 1925, proposed in her Ph.D. thesis an explanation for the composition of stars in terms of the relative abundances of hydrogen and helium. https://en.wikipedia.org/wiki/Cecilia_Payne-Gaposchkin
            "payne",

            // Radia Perlman is a software designer and network engineer and most famous for her invention of the spanning-tree protocol (STP). https://en.wikipedia.org/wiki/Radia_Perlman
            "perlman",

            // Thomas Gautier Pesquet - French astronaut. https://en.wikipedia.org/wiki/Thomas_Pesquet
            "pesquet",

            // Pablo Picasso was a Spanish painter, sculptor, printmaker, ceramicist, stage designer, poet and playwright. https://en.wikipedia.org/wiki/Pablo_Picasso
            "picasso",

            // Rob Pike was a key contributor to Unix, Plan 9, the X graphic system, utf-8, and the Go programming language. https://en.wikipedia.org/wiki/Rob_Pike
            "pike",

            // Joseph Plateau - Belgian physisist known for being one of the first persons to demonstrate the illusion of moving image. https://en.wikipedia.org/wiki/Joseph_Plateau
            "plateau",

            // Henri Poincaré made fundamental contributions in several fields of mathematics. https://en.wikipedia.org/wiki/Henri_Poincar%C3%A9
            "poincare",

            // Siméon Denis Poisson - French mathematician, geometer, and physicist, obtained many important results, was the final leading opponent of the wave theory of light and was proven wrong. https://en.wikipedia.org/wiki/Sim%C3%A9on_Denis_Poisson
            "poisson",

            // Laura Poitras is a director and producer whose work, made possible by open source crypto tools, advances the causes of truth and freedom of information by reporting disclosures by whistleblowers such as Edward Snowden. https://en.wikipedia.org/wiki/Laura_Poitras
            "poitras",

            // Claudius Ptolemy - a Greco-Egyptian writer of Alexandria, known as a mathematician, astronomer, geographer, astrologer, and poet of a single epigram in the Greek Anthology. https://en.wikipedia.org/wiki/Ptolemy
            "ptolemy",

            // C. V. Raman - Indian physicist who won the Nobel Prize in 1930 for proposing the Raman effect. https://en.wikipedia.org/wiki/C._V._Raman
            "raman",

            // Srinivasa Ramanujan - Indian mathematician and autodidact who made extraordinary contributions to mathematical analysis, number theory, infinite series, and continued fractions. https://en.wikipedia.org/wiki/Srinivasa_Ramanujan
            "ramanujan",

            // Sally Kristen Ride was an American physicist and astronaut. She was the first American woman in space, and the youngest American astronaut. https://en.wikipedia.org/wiki/Sally_Ride
            "ride",

            // Dennis Ritchie - co-creator of UNIX and the C programming language. https://en.wikipedia.org/wiki/Dennis_Ritchie
            "ritchie",

            // Wilhelm Conrad Röntgen - German physicist who was awarded the first Nobel Prize in Physics in 1901 for the discovery of X-rays (Röntgen rays). https://en.wikipedia.org/wiki/Wilhelm_R%C3%B6ntgen
            "roentgen",

            // Rosalind Franklin - English chemist and X-ray crystallographer whose research was critical to the understanding of DNA. https://en.wikipedia.org/wiki/Rosalind_Franklin
            "rosalind",

            // Peter Paul Rubens - Flemish baroque painter. https://en.wikipedia.org/wiki/Peter_Paul_Rubens
            "rubens",

            // Meghnad Saha - Indian astrophysicist best known for his development of the Saha equation, used to describe chemical and physical conditions in stars. https://en.wikipedia.org/wiki/Meghnad_Saha
            "saha",

            // Jean E. Sammet developed FORMAC, the first widely used computer language for symbolic manipulation of mathematical formulas. https://en.wikipedia.org/wiki/Jean_E._Sammet
            "sammet",

            // Frederick Sanger - British biochemist who worked on protein structure and DNA sequencing. https://en.wikipedia.org/wiki/Frederick_Sanger
            "sanger",

            // Adolphe Sax - inventor of the saxophone. https://en.wikipedia.org/wiki/Adolphe_Sax
            "sax",

            // Claude Shannon - The father of information theory and founder of digital circuit design theory. https://en.wikipedia.org/wiki/Claude_Shannon
            "shannon",

            // Carol Shaw - Originally an Atari employee, Carol Shaw is said to be the first female video game designer. https://en.wikipedia.org/wiki/Carol_Shaw_(video_game_designer)
            "shaw",

            // Dame Stephanie "Steve" Shirley - Founded a software company in 1962 employing women working from home. https://en.wikipedia.org/wiki/Steve_Shirley
            "shirley",

            // William Shockley co-invented the transistor. https://en.wikipedia.org/wiki/William_Shockley
            "shockley",

            // Françoise Barré-Sinoussi - French virologist and Nobel Prize Laureate in Physiology or Medicine; her work was fundamental in identifying HIV as the cause of AIDS. https://en.wikipedia.org/wiki/Fran%C3%A7oise_Barr%C3%A9-Sinoussi
            "sinoussi",

            // Betty Snyder - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Betty_Holberton
            "snyder",

            // Ernest Solvay - Belgian chemist. Founded Solvay & Cie and began a series of conferences in physics in 1911. https://en.wikipedia.org/wiki/Ernest_Solvay
            "solvay",

            // Frances Spence - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Frances_Spence
            "spence",

            // Richard Matthew Stallman - the founder of the Free Software movement, the GNU project, the Free Software Foundation, and the League for Programming Freedom. He also invented the concept of copyleft to protect the ideals of this movement, and enshrined this concept in the widely-used GPL (General Public License) for software. https://en.wikiquote.org/wiki/Richard_Stallman
            "stallman",

            // Arthur Harold Stone is credited with the discovery of the first flexagon. https://en.wikipedia.org/wiki/Arthur_Harold_Stone
            "stone",

            // Michael Stonebraker is a database research pioneer and architect of Ingres, Postgres, VoltDB and SciDB. Winner of 2014 ACM Turing Award. https://en.wikipedia.org/wiki/Michael_Stonebraker
            "stonebraker",

            // Janese Swanson (with others) developed the first of the Carmen Sandiego games. She went on to found Girl Tech. https://en.wikipedia.org/wiki/Janese_Swanson
            "swanson",

            // Aaron Swartz was influential in creating RSS, Markdown, Creative Commons, Reddit, and much of the internet as we know it today. He was devoted to freedom of information on the web. https://en.wikiquote.org/wiki/Aaron_Swartz
            "swartz",

            // Bertha Swirles was a theoretical physicist who made a number of contributions to early quantum theory. https://en.wikipedia.org/wiki/Bertha_Swirles
            "swirles",

            // Nikola Tesla invented the AC electric system and every gadget ever used by a James Bond villain. https://en.wikipedia.org/wiki/Nikola_Tesla
            "tesla",

            // Ken Thompson - co-creator of UNIX and the C programming language. https://en.wikipedia.org/wiki/Ken_Thompson
            "thompson",

            // Evangelista Torricelli was an Italian physicist and mathematician, best known for his invention of the barometer. https://en.wikipedia.org/wiki/Evangelista_Torricelli
            "torricelli",

            // Linus Torvalds invented Linux and Git. https://en.wikipedia.org/wiki/Linus_Torvalds
            "torvalds",

            // Bryant Tuckerman invented the Tuckerman traverse method for revealing all the faces of a flexagon. https://en.wikipedia.org/wiki/Bryant_Tuckerman
            "tuckerman",

            // Alan Turing was a founding father of computer science. https://en.wikipedia.org/wiki/Alan_Turing.
            "turing",

            // Varahamihira - Ancient Indian mathematician who discovered trigonometric formulae during 505-587 CE. https://en.wikipedia.org/wiki/Var%C4%81hamihira#Contributions
            "varahamihira",

            // Craig Venter - American biotechnologist, biochemist, geneticist, and businessman. Known for being involved with sequencing the second human genome. https://en.wikipedia.org/wiki/Craig_Venter
            "venter",

            // Sir Mokshagundam Visvesvaraya - is a notable Indian engineer.  He is a recipient of the Indian Republic's highest honour, the Bharat Ratna, in 1955. On his birthday, 15 September is celebrated as Engineer's Day in India in his memory. https://en.wikipedia.org/wiki/Visvesvaraya
            "visvesvaraya",

            // Christiane Nüsslein-Volhard - German biologist, won Nobel Prize in Physiology or Medicine in 1995 for research on the genetic control of embryonic development. https://en.wikipedia.org/wiki/Christiane_N%C3%BCsslein-Volhard
            "volhard",

            // Alessandro Giuseppe Antonio Anastasio Volta was an Italian physicist, chemist, and a pioneer of electricity and power. https://en.wikipedia.org/wiki/Alessandro_Volta
            "volta",

            // C. H. Waddington - British developmental biologist, paleontologist, geneticist, embryologist and philosopher who laid the foundations for systems biology, epigenetics, and evolutionary developmental biology. https://en.wikipedia.org/wiki/C._H._Waddington
            "waddington",

            // James Watson - American molecular biologist, geneticist and zoologist, best known as one of the co-discoverers of the structure of DNA. https://en.wikipedia.org/wiki/James_Watson
            "watson",

            // Marlyn Wescoff - one of the original programmers of the ENIAC. https://en.wikipedia.org/wiki/ENIAC - https://en.wikipedia.org/wiki/Marlyn_Meltzer
            "wescoff",

            // Andrew Wiles - Notable British mathematician who proved the enigmatic Fermat's Last Theorem. https://en.wikipedia.org/wiki/Andrew_Wiles
            "wiles",

            // Roberta Williams, did pioneering work in graphical adventure games for personal computers, particularly the King's Quest series. https://en.wikipedia.org/wiki/Roberta_Williams
            "williams",

            // Sophie Wilson designed the first Acorn Micro-Computer and the instruction set for ARM processors. https://en.wikipedia.org/wiki/Sophie_Wilson
            "wilson",

            // Jeannette Wing - co-developed the Liskov substitution principle. https://en.wikipedia.org/wiki/Jeannette_Wing
            "wing",

            // Carl Richard Woese - American microbiologist defined the Archaea kingdom of life by his pioneering phylogenetic taxonomy classification using 16S ribosomal RNA. https://en.wikipedia.org/wiki/Carl_Woese
            "woese",

            // Steve Wozniak invented the Apple I and Apple II. https://en.wikipedia.org/wiki/Steve_Wozniak
            "wozniak",

            // The Wright brothers, Orville and Wilbur - credited with inventing and building the world's first successful airplane and making the first controlled, powered and sustained heavier-than-air human flight. https://en.wikipedia.org/wiki/Wright_brothers
            "wright",

            // Rosalyn Sussman Yalow - Rosalyn Sussman Yalow was an American medical physicist, and a co-winner of the 1977 Nobel Prize in Physiology or Medicine for development of the radioimmunoassay technique. https://en.wikipedia.org/wiki/Rosalyn_Sussman_Yalow
            "yalow",

            // Ada Yonath - an Israeli crystallographer, the first woman from the Middle East to win a Nobel prize in the sciences. https://en.wikipedia.org/wiki/Ada_Yonath
            "yonath"

    ]

    private static Random RND = new Random()

    /**
    * @return A random generated name string
    */
    static String next() {
        def a = ADJECTIVES[ RND.nextInt(ADJECTIVES.size()) ]
        def b = NAMES[ RND.nextInt(NAMES.size()) ]
        return "${a}_${b}"
    }

    static String next(String... skip) {
        next(skip as List<String>)
    }

    static String next(Collection<String> skip) {
        while( true ) {
            final result = next()
            if( !skip.contains(result) )
                return result
        }
    }

}
