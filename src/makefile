# Название исполняемого файла
TARGET = a

# Компилятор и флаги
CXX = g++
CXXFLAGS = -w

# Список исходников
SOURCES = $(wildcard *.cc)

# Правило сборки
$(TARGET): $(SOURCES)
	$(CXX) -o $(TARGET) $(SOURCES) $(CXXFLAGS)

# Очистка
clean:
	rm -f $(TARGET)
